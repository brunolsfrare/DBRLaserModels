function [PumpAbsCrossSection, PumpEmCrossSection,LaserAbsCrossSection, LaserEmCrossSection,PumpBGLoss,LaserBGLoss,c,h,PumpFreq,LaserFreq,PumpWavelength, LaserWavelength,PumpPhotonEnergy,LaserPhotonEnergy, GratingLoss] = LoadParameters(dopant, PumpWavelength,LaserWavelength,PumpBGLossdB,LaserBGLossdB,GratingLossdB)
%Function that loads all simulation parameters to the main script

%Pump Absorption Cross Section
if dopant=='Er'
    if PumpWavelength==980
        %Not implemented
    else if PumpWavelength>=1460 && PumpWavelength<=1639
        %Loads Pump Abs Cross Section
        AbscsDat = 'TeO2_Abs_V1.dat';
        dat = importdata(AbscsDat);
        lambdaCS = dat(:,1); %(nm)
        sigma = dat(:,2); %(in cm^2)
            ind = find(lambdaCS > PumpWavelength, 1);
            lambdaA = lambdaCS(ind-1); lambdaB = lambdaCS(ind);
            sigmaA = sigma(ind-1); sigmaB = sigma(ind);
            PumpAbsCrossSection = sigmaA + (PumpWavelength - lambdaA) * ((sigmaB - sigmaA)/(lambdaB - lambdaA));

        %Loads Pump Em Cross Section
        EmcsDat = 'TeO2_EmASE_V1.dat';
        dat2 = importdata(EmcsDat);
        lambdaCS = dat2(:,1); %(nm)
        sigma = dat2(:,2); %(in cm^2)
            ind = find(lambdaCS > PumpWavelength, 1);
            lambdaA = lambdaCS(ind-1); lambdaB = lambdaCS(ind);
            sigmaA = sigma(ind-1); sigmaB = sigma(ind);
            PumpEmCrossSection = sigmaA + (PumpWavelength - lambdaA) * ((sigmaB - sigmaA)/(lambdaB - lambdaA));
    end
PumpAbsCrossSection = PumpAbsCrossSection*1E-4; %convert to m²
PumpEmCrossSection = PumpEmCrossSection*1E-4; %convert to m²
    end

        %Loads Laser Abs Cross Section
        AbscsDat = 'TeO2_Abs_V1.dat';
        dat = importdata(AbscsDat);
        lambdaCS = dat(:,1); %(nm)
        sigma = dat(:,2); %(in cm^2)
            ind = find(lambdaCS > LaserWavelength, 1);
            lambdaA = lambdaCS(ind-1); lambdaB = lambdaCS(ind);
            sigmaA = sigma(ind-1); sigmaB = sigma(ind);
            LaserAbsCrossSection = sigmaA + (LaserWavelength - lambdaA) * ((sigmaB - sigmaA)/(lambdaB - lambdaA));


        %Loads Laser Em Cross Section
        EmcsDat = 'TeO2_EmASE_V1.dat';
        dat2 = importdata(EmcsDat);
        lambdaCS = dat2(:,1); %(nm)
        sigma = dat2(:,2); %(in cm^2)
            ind = find(lambdaCS > LaserWavelength, 1);
            lambdaA = lambdaCS(ind-1); lambdaB = lambdaCS(ind);
            sigmaA = sigma(ind-1); sigmaB = sigma(ind);
            LaserEmCrossSection = sigmaA + (LaserWavelength - lambdaA) * ((sigmaB - sigmaA)/(lambdaB - lambdaA));


LaserAbsCrossSection = LaserAbsCrossSection*1E-4; %convert to m²
LaserEmCrossSection = LaserEmCrossSection*1E-4; %convert to m²          
end

%Convert losses from dB/cm to linear and meters (m^-1)
% "100*" converts to dB/m and 10*log10(exp(0)) converts to m^-1
PumpBGLoss= 100*PumpBGLossdB/(10*log10(exp(1))); 
LaserBGLoss=100*LaserBGLossdB/(10*log10(exp(1)));
GratingLoss = 100*GratingLossdB/(10*log10(exp(1)));

%Loads physical constants
c = physconst('LightSpeed'); %Speed of light in m/s
h = 6.62607015E-34; %Planck's constant in m²kg/s

%Convert wavelengths to m
PumpWavelength=PumpWavelength*1E-9;
LaserWavelength=LaserWavelength*1E-9;

%Pump and laser frequency calculation
PumpFreq = c/PumpWavelength;
LaserFreq = c/LaserWavelength;

%Photon energy
PumpPhotonEnergy = h*PumpFreq;
LaserPhotonEnergy = h*LaserFreq;
end