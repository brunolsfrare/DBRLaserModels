clear all;

%Mode areas and overlap with tellurite (from RSoft)
ModeArea_pump = 1.24551e-12;
ModeArea_laser = 1.61908e-12;
OverlapTellurite_pump = 0.633;
OverlapTellurite_laser = 0.573;

%Load simulation parameters (cross sections, rate equations, ...)
dopant = 'Tm'; %choose rare-earth: Er, Er-Yb, Tm, Yb,Pr,... Only Tm implemented in this version
PumpWavelength = 1610; %choose pump power. For Tm, a value between 1510 and 2100.
LaserWavelength = 1875; %in nm. For Tm, between 1565 and 2200

PumpBGLossdB = 0.5;%Losses(Lsweep);%0.5; %Pump background propagation loss in dB/cm
LaserBGLossdB = 0.5;% Losses(Lsweep);%0.5; %Laser wavelength background propagation loss in dB/cm
GratingLossdB = 0; %Grating loss in dB/cm 
%c is speed of light, h is Planck's constant
[PumpAbsCrossSection, PumpEmCrossSection,LaserAbsCrossSection, LaserEmCrossSection,PumpBGLoss,LaserBGLoss,c,h,PumpFreq,LaserFreq, PumpWavelength, LaserWavelength,PumpPhotonEnergy,LaserPhotonEnergy, GratingLoss] = LoadParameters(dopant,PumpWavelength,LaserWavelength,PumpBGLossdB,LaserBGLossdB,GratingLossdB);
tau1 = 0.3e-3; %Lifetime in seconds
tau3 = 345e-6; %Lifetime in seconds
BranchingRatio = 0.082;
W_ETU = 5E-18; %cm^3*s^-1. Will be converted to m³/s later
W_RETU = 1.8E-22; %in  m³/s
NT = 4E20; %Er concentration (cm-3). Will be converted to m-³ later
IncidentPumpPower_Fwd =[25, 50, 75, 100]*1E-3; %in W
IncidentPumpPower_Bwd = [0,0,0,0]*1E-3; %in W
Lin = 4E-3; %Input grating length in m
Lout = 0.5E-3; %Output grating length in m
Ltotal =2E-2; %Total device length (removing edge couplers) in m

W_ETU = W_ETU*1E-6; %Converts to m³/s
NT = NT*1E6; %Converts to m-³
kappa = 1150; %Grating coupling coefficient in m^-1
dZ = 0.01E-2; %Step size along z in m
NumOfSteps = int16(Ltotal/dZ+1); %Number of steps
Z = linspace(0,Ltotal,NumOfSteps); %converts steps to actual position along device

%Solve at the Bragg wavelength
dBeta = 0;
%Sets max iterations, correction, and tolerance
MaxIteration = 3E3;
PumpPowerError =ones(MaxIteration,1);
correction = 1E-7;
Tolerance = 1E-3;

%Initialize parameters
LaserPowerError=ones(numel(IncidentPumpPower_Fwd),NumOfSteps).*0.1;
ResultFwdLaserOutput = zeros(numel(IncidentPumpPower_Fwd),1);
ResultBwdLaserOutput = zeros(numel(IncidentPumpPower_Fwd),1);
AvgN0 = zeros(numel(IncidentPumpPower_Fwd),NumOfSteps);
AvgN1 = AvgN0;
AvgN3 = AvgN1;

%Loops that runs pump power sweep
for P=1:1:numel(IncidentPumpPower_Fwd)
    [FwdPumpPower,BwdPumpPower,FwdLaserPower,BwdLaserPower,N0,N1,N3,Gain_coeff,Abs_coeff,gain_coeff,abs_coeff,FwdLaserElectricField,BwdLaserElectricField] = CreateIterativeParameters(NumOfSteps,NT);
    FwdPumpPower(1)=IncidentPumpPower_Fwd(P);
    BwdPumpPower(1)=IncidentPumpPower_Bwd(P)/10;
    BwdPumpPower(NumOfSteps)=IncidentPumpPower_Bwd(P);
    FwdLaserPower(1) = 0;
    FwdLaserElectricField(1)=sqrt(FwdLaserPower(1));
    
    %Choose initial guess
    if P >1
        BwdLaserPower(1)=ResultBwdLaserOutput(P-1);
    else
        BwdLaserPower(1)=1E-6;
    end
        BwdLaserElectricField(1) = sqrt(BwdLaserPower(1));
    IterationN = 1;

    %main loop that propagates parameters
    while (abs(LaserPowerError(P,IterationN))>Tolerance && IterationN<MaxIteration) 

        BwdLaserPower(1)=BwdLaserPower(1)+correction;
        BwdLaserElectricField(1) = sqrt(BwdLaserPower(1));

        %Start guesses
        for z=1:1:NumOfSteps
            %Calculates parameters at z=0
            if z==1
                %Solve rate equations
                [N0(z),N1(z),N3(z)] = UpdatePouplations_Tm1610pump(tau1,tau3,W_ETU,FwdPumpPower(z)+BwdPumpPower(z),FwdLaserPower(z)+BwdLaserPower(z),NT,PumpAbsCrossSection, PumpEmCrossSection, PumpPhotonEnergy,LaserPhotonEnergy,LaserAbsCrossSection,LaserEmCrossSection,BranchingRatio,W_RETU,ModeArea_pump,ModeArea_laser,OverlapTellurite_pump,OverlapTellurite_laser);
                %Calculate gain and absorption coefficients
                Gain_coeff(z) = (LaserEmCrossSection.*N1(z)-LaserAbsCrossSection.*N0(z))*OverlapTellurite_laser;
                Abs_coeff(z) = (PumpEmCrossSection.*N1(z)-PumpAbsCrossSection.*N0(z))*OverlapTellurite_pump;
                
                BwdPumpPower(z) = BwdPumpPower(z+1)+(Abs_coeff(z)-PumpBGLoss)*BwdPumpPower(z+1)*dZ;

            elseif z>1 && Z(z)<=Lin %calculate parameters within input grating
                %Solve rate equations
                [N0(z),N1(z),N3(z)] = UpdatePouplations_Tm1610pump(tau1,tau3,W_ETU,FwdPumpPower(z-1)+BwdPumpPower(z-1),FwdLaserPower(z-1)+BwdLaserPower(z-1),NT,PumpAbsCrossSection, PumpEmCrossSection, PumpPhotonEnergy,LaserPhotonEnergy,LaserAbsCrossSection,LaserEmCrossSection,BranchingRatio,W_RETU,ModeArea_pump,ModeArea_laser,OverlapTellurite_pump,OverlapTellurite_laser);
                %Calculate gain and absorption coefficients
                Gain_coeff(z) = (LaserEmCrossSection.*N1(z-1)-LaserAbsCrossSection.*N0(z-1))*OverlapTellurite_laser;
                Abs_coeff(z) = (PumpEmCrossSection.*N1(z-1)-PumpAbsCrossSection.*N0(z-1))*OverlapTellurite_pump;

                %Propagates fields
                FwdLaserElectricField(z) = FwdLaserElectricField(z-1)+((-1i.*kappa*BwdLaserElectricField(z-1).*exp(+1i.*2.*dBeta.*dZ)) + (Gain_coeff(z-1)-LaserBGLoss-GratingLoss)*FwdLaserElectricField(z-1)).*dZ;
                BwdLaserElectricField(z) = BwdLaserElectricField(z-1)+((1i.*kappa*FwdLaserElectricField(z-1).*exp(-1i.*2.*dBeta.*dZ)) - (Gain_coeff(z-1)-LaserBGLoss-GratingLoss)*BwdLaserElectricField(z-1)).*dZ;
                
                %Calculates powers
                FwdLaserPower(z) = abs(FwdLaserElectricField(z))^2;
                BwdLaserPower(z) = abs(BwdLaserElectricField(z))^2;
                FwdPumpPower(z) = FwdPumpPower(z-1) +(Abs_coeff(z)-PumpBGLoss)*FwdPumpPower(z-1)*dZ;
                BwdPumpPower(z) = BwdPumpPower(z+1)+(Abs_coeff(z+1)-PumpBGLoss)*BwdPumpPower(z+1)*dZ;
               
            elseif Z(z)<=Ltotal-Lout %calculate parameters in region between gratings
                %Solve rate equations
                [N0(z),N1(z),N3(z)] = UpdatePouplations_Tm1610pump(tau1,tau3,W_ETU,FwdPumpPower(z-1)+BwdPumpPower(z-1),FwdLaserPower(z-1)+BwdLaserPower(z-1),NT,PumpAbsCrossSection, PumpEmCrossSection, PumpPhotonEnergy,LaserPhotonEnergy,LaserAbsCrossSection,LaserEmCrossSection,BranchingRatio,W_RETU,ModeArea_pump,ModeArea_laser,OverlapTellurite_pump,OverlapTellurite_laser);
                %Calculate gain and absorption coefficients
                Gain_coeff(z) = (LaserEmCrossSection.*N1(z-1)-LaserAbsCrossSection.*N0(z-1))*OverlapTellurite_laser;
                Abs_coeff(z) = (PumpEmCrossSection.*N1(z-1)-PumpAbsCrossSection.*N0(z-1))*OverlapTellurite_pump;

                %Propagates fields
                FwdLaserElectricField(z) = FwdLaserElectricField(z-1) + (Gain_coeff(z)-LaserBGLoss-GratingLoss)*FwdLaserElectricField(z-1).*dZ;
                BwdLaserElectricField(z) = BwdLaserElectricField(z-1) - (Gain_coeff(z)-LaserBGLoss-GratingLoss)*BwdLaserElectricField(z-1).*dZ;

                %Calculates powers
                FwdPumpPower(z) = FwdPumpPower(z-1) +(Abs_coeff(z)-PumpBGLoss)*FwdPumpPower(z-1)*dZ;
                BwdPumpPower(z) = BwdPumpPower(z+1)+(Abs_coeff(z)-PumpBGLoss)*BwdPumpPower(z+1)*dZ;
                FwdLaserPower(z) = abs(FwdLaserElectricField(z))^2;
                BwdLaserPower(z) = abs(BwdLaserElectricField(z))^2;

               
            elseif z<NumOfSteps %calculate parameters within output grating
                %Solve rate equations
                [N0(z),N1(z),N3(z)] = UpdatePouplations_Tm1610pump(tau1,tau3,W_ETU,FwdPumpPower(z-1)+BwdPumpPower(z-1),FwdLaserPower(z-1)+BwdLaserPower(z-1),NT,PumpAbsCrossSection, PumpEmCrossSection, PumpPhotonEnergy,LaserPhotonEnergy,LaserAbsCrossSection,LaserEmCrossSection,BranchingRatio,W_RETU,ModeArea_pump,ModeArea_laser,OverlapTellurite_pump,OverlapTellurite_laser);
                %Calculate gain and absorption coefficients for each mesh point
                Gain_coeff(z) = (LaserEmCrossSection.*N1(z-1)-LaserAbsCrossSection.*N0(z-1))*OverlapTellurite_laser;
                Abs_coeff(z) = (PumpEmCrossSection.*N1(z-1)-PumpAbsCrossSection.*N0(z-1))*OverlapTellurite_pump;
                
                %Propagates fields        
                FwdLaserElectricField(z) = FwdLaserElectricField(z-1)+((-1i.*kappa*BwdLaserElectricField(z-1).*exp(+1i.*2.*dBeta.*(-dZ))) + (Gain_coeff(z)-LaserBGLoss-GratingLoss)*FwdLaserElectricField(z-1)).*(-dZ);
                BwdLaserElectricField(z) = BwdLaserElectricField(z-1)+((1i.*kappa*FwdLaserElectricField(z-1).*exp(-1i.*2.*dBeta.*(-dZ))) - (Gain_coeff(z)-LaserBGLoss-GratingLoss)*BwdLaserElectricField(z-1)).*(-dZ);
    
                %Calculates powers
                FwdLaserPower(z) = abs(FwdLaserElectricField(z))^2;
                BwdLaserPower(z) = abs(BwdLaserElectricField(z))^2;
                FwdPumpPower(z) = FwdPumpPower(z-1) +(Abs_coeff(z)-PumpBGLoss)*FwdPumpPower(z-1)*dZ;
                BwdPumpPower(z) = BwdPumpPower(z+1)+(Abs_coeff(z)-PumpBGLoss)*BwdPumpPower(z+1)*dZ;
              

            else
                %Solve rate equations
                [N0(z),N1(z),N3(z)] = UpdatePouplations_Tm1610pump(tau1,tau3,W_ETU,FwdPumpPower(z-1)+BwdPumpPower(z-1),FwdLaserPower(z-1)+BwdLaserPower(z-1),NT,PumpAbsCrossSection, PumpEmCrossSection, PumpPhotonEnergy,LaserPhotonEnergy,LaserAbsCrossSection,LaserEmCrossSection,BranchingRatio,W_RETU,ModeArea_pump,ModeArea_laser,OverlapTellurite_pump,OverlapTellurite_laser);
                %Calculate gain and absorption coefficients for each mesh point
                Gain_coeff(z) = (LaserEmCrossSection.*N1(z-1)-LaserAbsCrossSection.*N0(z-1))*OverlapTellurite_laser;
                Abs_coeff(z) = (PumpEmCrossSection.*N1(z-1)-PumpAbsCrossSection.*N0(z-1))*OverlapTellurite_pump;
                
                %Propagates fields
                FwdLaserElectricField(z) = FwdLaserElectricField(z-1)+((-1i.*kappa*BwdLaserElectricField(z-1).*exp(+1i.*2.*dBeta.*dZ)) + (Gain_coeff(z)-LaserBGLoss-GratingLoss)*FwdLaserElectricField(z-1)).*dZ;
                BwdLaserElectricField(z) = BwdLaserElectricField(z-1)+((1i.*kappa*FwdLaserElectricField(z-1).*exp(-1i.*2.*dBeta.*dZ)) - (Gain_coeff(z)-LaserBGLoss-GratingLoss)*BwdLaserElectricField(z-1)).*dZ;
    
                %Calculates powers
                FwdLaserPower(z) = abs(FwdLaserElectricField(z))^2;
                BwdLaserPower(z) = abs(BwdLaserElectricField(z))^2;
                FwdPumpPower(z) = FwdPumpPower(z-1) +(Abs_coeff(z)-PumpBGLoss)*FwdPumpPower(z-1)*dZ;
            end

        end
    %Check boundary conditions and update values
    IterationN=IterationN+1;
    LaserPowerError(P,IterationN) =BwdLaserElectricField(NumOfSteps);
   
    disp((BwdLaserElectricField(NumOfSteps)));

    if IterationN==MaxIteration
        disp("MaxIteration reached, no convergence");
        FwdLaserPower(:)=0;
        BwdLaserPower(:)=0;
    end

    end
    ResultFwdLaserOutput(P)=FwdLaserPower(NumOfSteps);
    ResultBwdLaserOutput(P)=BwdLaserPower(1);

    fprintf("Finished %d",1000*IncidentPumpPower_Fwd(P));
end
a=linspace(1,IterationN,IterationN);

filename = sprintf("Test");
save(filename,"ResultFwdLaserOutput","ResultBwdLaserOutput","IncidentPumpPower_Fwd","IncidentPumpPower_Bwd","AvgN0","AvgN1","AvgN3","Z","FwdLaserPower","BwdLaserPower","FwdPumpPower","BwdPumpPower","Gain_coeff","Abs_coeff","FwdLaserElectricField","BwdLaserElectricField");

