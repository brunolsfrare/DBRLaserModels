function [FwdPumpPower,BwdPumpPower,FwdLaserPower,BwdLaserPower,N0,N1,N2,N0_active,N1_active,N2_active,N0_quench,N1_quench,N2_quench,Gain_coeff,Abs_coeff,gain_coeff,abs_coeff,FwdLaserElectricField,BwdLaserElectricField] = CreateIterativeParameters(RefractiveIndex,NumOfSteps,GainMediumRefractiveIndex,NT,NActive,NQuench)

%   Creates the variables that are updated at each step of the simulation

%Create total, active, quenched populations, gain and abs coefficients
N0 = zeros(NumOfSteps,numel(RefractiveIndex(:,1)),numel(RefractiveIndex(1,:)));
N1 = zeros(NumOfSteps,numel(RefractiveIndex(:,1)),numel(RefractiveIndex(1,:)));
N2 = zeros(NumOfSteps,numel(RefractiveIndex(:,1)),numel(RefractiveIndex(1,:)));
N0_active = zeros(NumOfSteps,numel(RefractiveIndex(:,1)),numel(RefractiveIndex(1,:)));
N1_active = zeros(NumOfSteps,numel(RefractiveIndex(:,1)),numel(RefractiveIndex(1,:)));
N2_active = zeros(NumOfSteps,numel(RefractiveIndex(:,1)),numel(RefractiveIndex(1,:)));
N0_quench = zeros(NumOfSteps,numel(RefractiveIndex(:,1)),numel(RefractiveIndex(1,:)));
N1_quench = zeros(NumOfSteps,numel(RefractiveIndex(:,1)),numel(RefractiveIndex(1,:)));
N2_quench = zeros(NumOfSteps,numel(RefractiveIndex(:,1)),numel(RefractiveIndex(1,:)));
gain_coeff = zeros(NumOfSteps,numel(RefractiveIndex(:,1)),numel(RefractiveIndex(1,:)));
abs_coeff = zeros(NumOfSteps,numel(RefractiveIndex(:,1)),numel(RefractiveIndex(1,:)));

%Initialize powers/field and average coeffieicnets as zero
FwdPumpPower = zeros(NumOfSteps,1).*1E-8; 
BwdPumpPower = zeros(NumOfSteps,1).*1E-8;
FwdLaserPower = zeros(NumOfSteps,1).*1E-8;
BwdLaserPower = zeros(NumOfSteps,1).*1E-8;
Gain_coeff = zeros(NumOfSteps,1);
Abs_coeff = zeros(NumOfSteps,1);
FwdLaserElectricField = FwdLaserPower;
BwdLaserElectricField = BwdLaserPower;

%Initialize total, active, and quenched ion populations in the ground state
%in the mesh points that correspond to TeO2
for i=1:1:numel(N0(1,:,1))
    for j=1:1:numel(N0(1,1,:))
        if RefractiveIndex(i,j) == GainMediumRefractiveIndex
            N0(:,i,j) = NT;
            N0_active(:,i,j) = NActive;
            N0_quench(:,i,j) = NQuench;

        end
    end
end
end