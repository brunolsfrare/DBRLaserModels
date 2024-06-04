function [Areas,RefractiveIndex,PumpNormIdist,LaserNormIdist] = GenerateMatrices(mesh,RefractiveIndex,PumpElectricFieldProfile,LaserElectricFieldProfile)
%Function that receives RSoft files and generate matrices for simulations


RefractiveIndex(:,1) = [];
LaserElectricFieldProfile(1,:) = [];
PumpElectricFieldProfile(1,:) = [];


%Mesh element area (in this case, 10 x 10 nm uniform mesh)
Areas=1E-16;

%Create intensity profiles from electric field profiles
PumpIntensityProfile=PumpElectricFieldProfile.^2;
LaserIntensityProfile=LaserElectricFieldProfile.^2;

%Create normalized intensity distributions
PumpNormIdist = (PumpIntensityProfile.*Areas)./(sum(sum(PumpIntensityProfile.*Areas)));
LaserNormIdist = (LaserIntensityProfile.*Areas)./(sum(sum(LaserIntensityProfile.*Areas)));

%Arrange distributions in the right orientation
RefractiveIndex=transpose(RefractiveIndex);
LaserNormIdist = transpose(LaserNormIdist);
PumpNormIdist = transpose(PumpNormIdist);

end