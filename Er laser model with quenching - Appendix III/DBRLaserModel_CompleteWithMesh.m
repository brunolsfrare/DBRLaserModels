clear all;
mesh=readtable("T92Mesh.txt"); %convert to txt before importing (MSH file)
RefractiveIndex = importdata("T92RefractiveIndex.txt"); %convert to txt and delete header before importing
PumpElectricFieldProfile = importdata("T92PumpProfile.txt"); %In RSoft, need to go to Mode solver>output>change to amplitude. It has to be a uniform mesh or RSoft will not save it with the same number of points as the other files (Ex.m00). Don't make sim window too large to make the sim a bit faster
LaserElectricFieldProfile = importdata("T92LaserProfile.txt"); %it has to be a uniform mesh or RSoft will not save it with the same number of points as the other files (Ex.m00)
[Areas,RefractiveIndex,PumpNormIdist,LaserNormIdist] = GenerateMatrices(mesh,RefractiveIndex,PumpElectricFieldProfile,LaserElectricFieldProfile);
LaserNormEdist = sqrt(LaserNormIdist)./sum(sum(sqrt(LaserNormIdist)));
%CHECK ALL UNITS: cross sections, wavelengths....make sure everything is in
%meters,W, and SI

%Load simulation parameters (cross sections, rate equations, ...)
dopant = 'Er'; %choose rare-earth: Er, Er-Yb, Tm, Yb,Pr,... Only Er implemented so far
PumpWavelength = 1470; %choose pump wavelength. For Er, a value between 1460 and 1639. 980 has not been implemented.
LaserWavelength = 1550; %in nm. For Er, between 1460 and 1640
 
PumpBGLossdB = 0.5; %Pump background propagation loss in dB/cm
LaserBGLossdB = 0.5; %Laser wavelength background propagation loss in dB/cm
GratingLossdB = 0; %Grating loss in dB/cm 
%c is speed of light, h is Planck's constant
[PumpAbsCrossSection, PumpEmCrossSection,LaserAbsCrossSection, LaserEmCrossSection,PumpBGLoss,LaserBGLoss,c,h,PumpFreq,LaserFreq, PumpWavelength, LaserWavelength,PumpPhotonEnergy,LaserPhotonEnergy, GratingLoss] = LoadParameters(dopant,PumpWavelength,LaserWavelength,PumpBGLossdB,LaserBGLossdB,GratingLossdB);
tau1 = 0.62e-3; %Excited-state lifetime of active ion in seconds
tau1_QuenchedIons = 1e-6; %Excited-state lifetime of quenched ions in seconds
tau2 = 1e-6; %Lifetime of energy level E2 in seconds
W_ETU = 2.8E-18; %cm^3*s^-1, from Henry Frankis' thesis. Will be converted to m³/s later
NT = 2.5E20; %Er concentration (cm-3). Will be converted to m-³ later
IncidentPumpPower_Fwd =[25, 50, 75, 100]*1E-3; %in W
IncidentPumpPower_Bwd = [0,0,0,0]*1E-3; %in W
Lin = 4E-3; %Input grating length in m
Lout = 0.85E-3; %Output grating length in m
Ltotal = 2E-2; %Total device length (removing edge couplers) in m

GainMediumRefractiveIndex = 2.05; %MUST match value used in RSoft (RefractiveIndex matrix)
W_ETU = W_ETU*1E-6; %Converts to m³/s
NT = NT*1E6; %Converts to m-³
Quenching_fraction = 0.2;
NQuench = NT*Quenching_fraction;
NActive = NT-NQuench;
kappa = 1150; %Grating coupling coefficient in m^-1
dZ = 0.01E-2; %Stepsize of our sims along Z in m
NumOfSteps = int16(Ltotal/dZ+1); %Number of steps
Z = linspace(0,Ltotal,NumOfSteps); %converts steps to actual position along device

%sets simulation to Bragg wavelength
dBeta = 0;
%Sets tolerance, max iterations, guess correction
MaxIteration = 3E3;
PumpPowerError =ones(MaxIteration,1);
correction = 1E-7; %usually 1E-8 or 1E-7
Tolerance = 1E-3;

%Initializes variables
LaserPowerError=ones(numel(IncidentPumpPower_Fwd),NumOfSteps).*0.1;
ResultFwdLaserOutput = zeros(numel(IncidentPumpPower_Fwd),1);
ResultBwdLaserOutput = zeros(numel(IncidentPumpPower_Fwd),1);
AvgN0 = zeros(numel(IncidentPumpPower_Fwd),NumOfSteps);
AvgN1 = AvgN0;
AvgN2 = AvgN1;

%Loop that runs pump power sweep
for P=1:1:numel(IncidentPumpPower_Fwd)
    [FwdPumpPower,BwdPumpPower,FwdLaserPower,BwdLaserPower,N0,N1,N2,N0_active,N1_active,N2_active,N0_quench,N1_quench,N2_quench,Gain_coeff,Abs_coeff,gain_coeff,abs_coeff,FwdLaserElectricField,BwdLaserElectricField] = CreateIterativeParameters(RefractiveIndex,NumOfSteps,GainMediumRefractiveIndex,NT,NActive,NQuench);
    %Initialize powers
    FwdPumpPower(1)=IncidentPumpPower_Fwd(P);
    BwdPumpPower(1)=IncidentPumpPower_Bwd(P)/10;
    BwdPumpPower(NumOfSteps)=IncidentPumpPower_Bwd(P);
    FwdLaserPower(1) = 0;
    FwdLaserElectricField(1)=sqrt(FwdLaserPower(1));
    
    %If previous pump power converged, use the result as initial guess
    if P >1 && ResultBwdLaserOutput(P-1)~=0
        BwdLaserPower(1)=ResultBwdLaserOutput(P-1);
    else
        %otherwise, choose initial guess (usually between 1E-8 and 5E-6)
        BwdLaserPower(1)=1E-7;
    end
        BwdLaserElectricField(1) = sqrt(BwdLaserPower(1));
    IterationN = 1;

    %Main loop that propagates all parameters along Z
    while (abs(LaserPowerError(P,IterationN))>Tolerance && IterationN<MaxIteration) 

        %Update initial guess
        BwdLaserPower(1)=BwdLaserPower(1)+correction;
        BwdLaserElectricField(1) = sqrt(BwdLaserPower(1));

        
        %Start guesses
        for z=1:1:NumOfSteps
            %Calculate paramters at first point
            if z==1
                %Update populations
                [N0_active(z,:,:),N1_active(z,:,:),N2_active(z,:,:)] = UpdatePouplations_Er1480pump(tau1,tau2,W_ETU,FwdPumpPower(z)+BwdPumpPower(z),FwdLaserPower(z)+BwdLaserPower(z),NActive,PumpNormIdist,PumpAbsCrossSection, PumpEmCrossSection, PumpPhotonEnergy,LaserPhotonEnergy,LaserNormIdist,LaserAbsCrossSection,LaserEmCrossSection,RefractiveIndex,GainMediumRefractiveIndex,Areas);
                [N0_quench(z,:,:),N1_quench(z,:,:),N2_quench(z,:,:)] = UpdatePouplations_Er1480pump(tau1_QuenchedIons,tau2,W_ETU,FwdPumpPower(z)+BwdPumpPower(z),FwdLaserPower(z)+BwdLaserPower(z),NQuench,PumpNormIdist,PumpAbsCrossSection, PumpEmCrossSection, PumpPhotonEnergy,LaserPhotonEnergy,LaserNormIdist,LaserAbsCrossSection,LaserEmCrossSection,RefractiveIndex,GainMediumRefractiveIndex,Areas);
                N0(z,:,:)=N0_active(z,:,:)+N0_quench(z,:,:);
                N1(z,:,:)=N1_active(z,:,:)+N1_quench(z,:,:);
                N2(z,:,:)=N2_active(z,:,:)+N2_quench(z,:,:);
                %Calculate gain and absorption coefficients for each mesh point
                gain_coeff(z,:,:) = LaserEmCrossSection.*N1(z,:,:)-LaserAbsCrossSection.*N0(z,:,:);
                abs_coeff(z,:,:) = PumpEmCrossSection.*N1(z,:,:)-PumpAbsCrossSection.*N0(z,:,:);
                %Calculate total (average) gain an absorption coefficients
                Gain_coeff(z) = sum(sum(LaserNormIdist.*squeeze(gain_coeff(z,:,:)))); 
                Abs_coeff(z) = sum(sum(PumpNormIdist.*squeeze(abs_coeff(z,:,:)))); 

                BwdPumpPower(z) = BwdPumpPower(z+1)+(Abs_coeff(z)-PumpBGLoss)*BwdPumpPower(z+1)*dZ;


            elseif z>1 && Z(z)<=Lin %calculate propagation within input grating
                %Update populations
                [N0_active(z,:,:),N1_active(z,:,:),N2_active(z,:,:)] = UpdatePouplations_Er1480pump(tau1,tau2,W_ETU,FwdPumpPower(z-1)+BwdPumpPower(z-1),FwdLaserPower(z-1)+BwdLaserPower(z-1),NActive,PumpNormIdist,PumpAbsCrossSection, PumpEmCrossSection, PumpPhotonEnergy,LaserPhotonEnergy,LaserNormIdist,LaserAbsCrossSection,LaserEmCrossSection,RefractiveIndex,GainMediumRefractiveIndex,Areas);
                [N0_quench(z,:,:),N1_quench(z,:,:),N2_quench(z,:,:)] = UpdatePouplations_Er1480pump(tau1_QuenchedIons,tau2,W_ETU,FwdPumpPower(z-1)+BwdPumpPower(z-1),FwdLaserPower(z-1)+BwdLaserPower(z-1),NQuench,PumpNormIdist,PumpAbsCrossSection, PumpEmCrossSection, PumpPhotonEnergy,LaserPhotonEnergy,LaserNormIdist,LaserAbsCrossSection,LaserEmCrossSection,RefractiveIndex,GainMediumRefractiveIndex,Areas);
                N0(z,:,:)=N0_active(z,:,:)+N0_quench(z,:,:);
                N1(z,:,:)=N1_active(z,:,:)+N1_quench(z,:,:);
                N2(z,:,:)=N2_active(z,:,:)+N2_quench(z,:,:);
                %Calculate gain and absorption coefficients for each mesh point
                gain_coeff(z,:,:) = LaserEmCrossSection.*N1(z-1,:,:)-LaserAbsCrossSection.*N0(z-1,:,:);
                abs_coeff(z,:,:) = PumpEmCrossSection.*N1(z-1,:,:)-PumpAbsCrossSection.*N0(z-1,:,:);
                %Calculate total gain an absorption coefficients
                Gain_coeff(z) = sum(sum(LaserNormIdist.*squeeze(gain_coeff(z,:,:)))); 
                Abs_coeff(z) =sum(sum(PumpNormIdist.*squeeze(abs_coeff(z,:,:))));
                
                %Propagates fields
                FwdLaserElectricField(z) = FwdLaserElectricField(z-1)+((-1i.*kappa*BwdLaserElectricField(z-1).*exp(+1i.*2.*dBeta.*dZ)) + (Gain_coeff(z)-LaserBGLoss-GratingLoss)*FwdLaserElectricField(z-1)).*dZ;
                BwdLaserElectricField(z) = BwdLaserElectricField(z-1)+((1i.*kappa*FwdLaserElectricField(z-1).*exp(-1i.*2.*dBeta.*dZ)) - (Gain_coeff(z)-LaserBGLoss-GratingLoss)*BwdLaserElectricField(z-1)).*dZ;
                
                %Calculates powers
                FwdLaserPower(z) = abs(FwdLaserElectricField(z))^2;
                BwdLaserPower(z) = abs(BwdLaserElectricField(z))^2;
                FwdPumpPower(z) = FwdPumpPower(z-1) +(Abs_coeff(z)-PumpBGLoss)*FwdPumpPower(z-1)*dZ;
                BwdPumpPower(z) = BwdPumpPower(z+1)+(Abs_coeff(z)-PumpBGLoss)*BwdPumpPower(z+1)*dZ;
               
            %Propagates parameters in region between gratings
            elseif Z(z)<=Ltotal-Lout
                %Update populations
                [N0_active(z,:,:),N1_active(z,:,:),N2_active(z,:,:)] = UpdatePouplations_Er1480pump(tau1,tau2,W_ETU,FwdPumpPower(z-1)+BwdPumpPower(z-1),FwdLaserPower(z-1)+BwdLaserPower(z-1),NActive,PumpNormIdist,PumpAbsCrossSection, PumpEmCrossSection, PumpPhotonEnergy,LaserPhotonEnergy,LaserNormIdist,LaserAbsCrossSection,LaserEmCrossSection,RefractiveIndex,GainMediumRefractiveIndex,Areas);
                [N0_quench(z,:,:),N1_quench(z,:,:),N2_quench(z,:,:)] = UpdatePouplations_Er1480pump(tau1_QuenchedIons,tau2,W_ETU,FwdPumpPower(z-1)+BwdPumpPower(z-1),FwdLaserPower(z-1)+BwdLaserPower(z-1),NQuench,PumpNormIdist,PumpAbsCrossSection, PumpEmCrossSection, PumpPhotonEnergy,LaserPhotonEnergy,LaserNormIdist,LaserAbsCrossSection,LaserEmCrossSection,RefractiveIndex,GainMediumRefractiveIndex,Areas);
                N0(z,:,:)=N0_active(z,:,:)+N0_quench(z,:,:);
                N1(z,:,:)=N1_active(z,:,:)+N1_quench(z,:,:);
                N2(z,:,:)=N2_active(z,:,:)+N2_quench(z,:,:);
                %Calculate gain and absorption coefficients for each mesh point
                gain_coeff(z,:,:) = LaserEmCrossSection.*N1(z-1,:,:)-LaserAbsCrossSection.*N0(z-1,:,:);
                abs_coeff(z,:,:) = PumpEmCrossSection.*N1(z-1,:,:)-PumpAbsCrossSection.*N0(z-1,:,:);
                %Calculate total gain an absorption coefficients
                Gain_coeff(z) = sum(sum(LaserNormIdist.*squeeze(gain_coeff(z,:,:)))); 
                Abs_coeff(z) = sum(sum(PumpNormIdist.*squeeze(abs_coeff(z,:,:)))); 
                
                %Propagates fields
                FwdLaserElectricField(z) = FwdLaserElectricField(z-1) + (Gain_coeff(z)-LaserBGLoss-GratingLoss)*FwdLaserElectricField(z-1).*dZ;
                BwdLaserElectricField(z) = BwdLaserElectricField(z-1) - (Gain_coeff(z)-LaserBGLoss-GratingLoss)*BwdLaserElectricField(z-1).*dZ;
                
                %Calculate powers
                FwdPumpPower(z) = FwdPumpPower(z-1) +(Abs_coeff(z)-PumpBGLoss)*FwdPumpPower(z-1)*dZ;
                BwdPumpPower(z) = BwdPumpPower(z+1)+(Abs_coeff(z)-PumpBGLoss)*BwdPumpPower(z+1)*dZ;
                FwdLaserPower(z) = abs(FwdLaserElectricField(z))^2;
                BwdLaserPower(z) = abs(BwdLaserElectricField(z))^2;

               %Propagate parameters within output grating region
            elseif z<NumOfSteps
                [N0_active(z,:,:),N1_active(z,:,:),N2_active(z,:,:)] = UpdatePouplations_Er1480pump(tau1,tau2,W_ETU,FwdPumpPower(z-1)+BwdPumpPower(z-1),FwdLaserPower(z-1)+BwdLaserPower(z-1),NActive,PumpNormIdist,PumpAbsCrossSection, PumpEmCrossSection, PumpPhotonEnergy,LaserPhotonEnergy,LaserNormIdist,LaserAbsCrossSection,LaserEmCrossSection,RefractiveIndex,GainMediumRefractiveIndex,Areas);
                [N0_quench(z,:,:),N1_quench(z,:,:),N2_quench(z,:,:)] = UpdatePouplations_Er1480pump(tau1_QuenchedIons,tau2,W_ETU,FwdPumpPower(z-1)+BwdPumpPower(z-1),FwdLaserPower(z-1)+BwdLaserPower(z-1),NQuench,PumpNormIdist,PumpAbsCrossSection, PumpEmCrossSection, PumpPhotonEnergy,LaserPhotonEnergy,LaserNormIdist,LaserAbsCrossSection,LaserEmCrossSection,RefractiveIndex,GainMediumRefractiveIndex,Areas);
                N0(z,:,:)=N0_active(z,:,:)+N0_quench(z,:,:);
                N1(z,:,:)=N1_active(z,:,:)+N1_quench(z,:,:);
                N2(z,:,:)=N2_active(z,:,:)+N2_quench(z,:,:);
                %Calculate gain and absorption coefficients for each mesh point
                gain_coeff(z,:,:) = LaserEmCrossSection.*N1(z-1,:,:)-LaserAbsCrossSection.*N0(z-1,:,:);
                abs_coeff(z,:,:) = PumpEmCrossSection.*N1(z-1,:,:)-PumpAbsCrossSection.*N0(z-1,:,:);
                %Calculate total gain an absorption coefficients
                Gain_coeff(z) = sum(sum(LaserNormIdist.*squeeze(gain_coeff(z,:,:))));
                Abs_coeff(z) = sum(sum(PumpNormIdist.*squeeze(abs_coeff(z,:,:)))); 
                
                %Propagates fields          
                FwdLaserElectricField(z) = FwdLaserElectricField(z-1)+((-1i.*kappa*BwdLaserElectricField(z-1).*exp(+1i.*2.*dBeta.*(-dZ))) + (Gain_coeff(z)-LaserBGLoss-GratingLoss)*FwdLaserElectricField(z-1)).*(-dZ);
                BwdLaserElectricField(z) = BwdLaserElectricField(z-1)+((1i.*kappa*FwdLaserElectricField(z-1).*exp(-1i.*2.*dBeta.*(-dZ))) - (Gain_coeff(z)-LaserBGLoss-GratingLoss)*BwdLaserElectricField(z-1)).*(-dZ);
                
                %Calculates powers
                FwdLaserPower(z) = abs(FwdLaserElectricField(z))^2;
                BwdLaserPower(z) = abs(BwdLaserElectricField(z))^2;
                FwdPumpPower(z) = FwdPumpPower(z-1) +(Abs_coeff(z)-PumpBGLoss)*FwdPumpPower(z-1)*dZ;
                BwdPumpPower(z) = BwdPumpPower(z+1)+(Abs_coeff(z)-PumpBGLoss)*BwdPumpPower(z+1)*dZ;
              
            %Calculate parameters in the last mesh point along Z
            else
                [N0_active(z,:,:),N1_active(z,:,:),N2_active(z,:,:)] = UpdatePouplations_Er1480pump(tau1,tau2,W_ETU,FwdPumpPower(z-1)+BwdPumpPower(z-1),FwdLaserPower(z-1)+BwdLaserPower(z-1),NActive,PumpNormIdist,PumpAbsCrossSection, PumpEmCrossSection, PumpPhotonEnergy,LaserPhotonEnergy,LaserNormIdist,LaserAbsCrossSection,LaserEmCrossSection,RefractiveIndex,GainMediumRefractiveIndex,Areas);
                [N0_quench(z,:,:),N1_quench(z,:,:),N2_quench(z,:,:)] = UpdatePouplations_Er1480pump(tau1_QuenchedIons,tau2,W_ETU,FwdPumpPower(z-1)+BwdPumpPower(z-1),FwdLaserPower(z-1)+BwdLaserPower(z-1),NQuench,PumpNormIdist,PumpAbsCrossSection, PumpEmCrossSection, PumpPhotonEnergy,LaserPhotonEnergy,LaserNormIdist,LaserAbsCrossSection,LaserEmCrossSection,RefractiveIndex,GainMediumRefractiveIndex,Areas);
                N0(z,:,:)=N0_active(z,:,:)+N0_quench(z,:,:);
                N1(z,:,:)=N1_active(z,:,:)+N1_quench(z,:,:);
                N2(z,:,:)=N2_active(z,:,:)+N2_quench(z,:,:);

                %Calculate gain and absorption coefficients for each mesh point
                gain_coeff(z,:,:) = LaserEmCrossSection.*N1(z-1,:,:)-LaserAbsCrossSection.*N0(z-1,:,:);
                abs_coeff(z,:,:) = PumpEmCrossSection.*N1(z-1,:,:)-PumpAbsCrossSection.*N0(z-1,:,:);
                %Calculate total gain an absorption coefficients
                Gain_coeff(z) = sum(sum(LaserNormIdist.*squeeze(gain_coeff(z,:,:)))); 
                Abs_coeff(z) = sum(sum(PumpNormIdist.*squeeze(abs_coeff(z,:,:)))); 
               
                %Propagates fields
                FwdLaserElectricField(z) = FwdLaserElectricField(z-1)+((-1i.*kappa*BwdLaserElectricField(z-1).*exp(+1i.*2.*dBeta.*dZ)) + (Gain_coeff(z)-LaserBGLoss-GratingLoss)*FwdLaserElectricField(z-1)).*dZ;
                BwdLaserElectricField(z) = BwdLaserElectricField(z-1)+((1i.*kappa*FwdLaserElectricField(z-1).*exp(-1i.*2.*dBeta.*dZ)) - (Gain_coeff(z)-LaserBGLoss-GratingLoss)*BwdLaserElectricField(z-1)).*dZ;

                %Calculates powers
                FwdLaserPower(z) = abs(FwdLaserElectricField(z))^2;
                BwdLaserPower(z) = abs(BwdLaserElectricField(z))^2;
                FwdPumpPower(z) = FwdPumpPower(z-1) +(Abs_coeff(z)-PumpBGLoss)*FwdPumpPower(z-1)*dZ;
            end
        %Calculates average populations along Z
        AvgN0(P,z) = mean(nonzeros(N0(z,:,:)));
        AvgN1(P,z) = mean(nonzeros(N1(z,:,:)));
        AvgN2(P,z) = mean(nonzeros(N2(z,:,:)));
        end
    %Check boundary conditions and update values
    IterationN=IterationN+1;
    LaserPowerError(P,IterationN) =BwdLaserElectricField(NumOfSteps);
   
    %Displays the value of Bwd Laser Electric field at the last point after
    %each iteration. This value should be negative initially and then as it
    %reaches zero, the code will converge. If it shows positive values,
    %choose smaller guess and correction. If it's negative but it's taking 
    %too long to converge, increase guess and correction values.If no matter
    %how low you choose your guess it always converges in the first 
    %iteratction, the device does not lase with that pump power
    disp((BwdLaserElectricField(NumOfSteps)));
    
    %If the model did not converge, no lasing
    if IterationN==MaxIteration
        disp("MaxIteration reached, no convergence");
        FwdLaserPower(:)=0;
        BwdLaserPower(:)=0;
    end

    end
    %Laser output powers after converging
    ResultFwdLaserOutput(P)=FwdLaserPower(NumOfSteps);
    ResultBwdLaserOutput(P)=BwdLaserPower(1);

    fprintf("Finished %d",1000*IncidentPumpPower_Fwd(P));
end
a=linspace(1,IterationN,IterationN);

filename = sprintf("Test");
save(filename,"ResultFwdLaserOutput","ResultBwdLaserOutput","IncidentPumpPower_Fwd","IncidentPumpPower_Bwd","AvgN0","AvgN1","AvgN2","Z","FwdLaserPower","BwdLaserPower","FwdPumpPower","BwdPumpPower","Gain_coeff","Abs_coeff","FwdLaserElectricField","BwdLaserElectricField","LaserPowerError");
                        
