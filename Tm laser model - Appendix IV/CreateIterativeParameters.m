function [FwdPumpPower,BwdPumpPower,FwdLaserPower,BwdLaserPower,N0,N1,N2,Gain_coeff,Abs_coeff,gain_coeff,abs_coeff,FwdLaserElectricField,BwdLaserElectricField] = CreateIterativeParameters(NumOfSteps,NT)

%   Creates the variables that are updated at each step of the simulation

N0 = ones(NumOfSteps,1)*NT;
N1 = zeros(NumOfSteps,1);
N2 = zeros(NumOfSteps,1);
gain_coeff = zeros(NumOfSteps,1);
abs_coeff = zeros(NumOfSteps,1);

FwdPumpPower = zeros(NumOfSteps,1).*1E-8;
BwdPumpPower = zeros(NumOfSteps,1).*1E-8;
FwdLaserPower = zeros(NumOfSteps,1).*1E-8;
BwdLaserPower = zeros(NumOfSteps,1).*1E-8;
Gain_coeff = zeros(NumOfSteps,1);
Abs_coeff = zeros(NumOfSteps,1);
FwdLaserElectricField = FwdLaserPower;
BwdLaserElectricField = BwdLaserPower;
