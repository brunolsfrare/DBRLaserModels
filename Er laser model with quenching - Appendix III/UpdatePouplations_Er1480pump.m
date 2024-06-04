function [N0,N1,N2] = UpdatePouplations_Er1480pump(tau1,tau2,W_ETU,PumpPower,LaserPower,NT,PumpNormIdist,PumpAbsCrossSection, PumpEmCrossSection, PumpPhotonEnergy,LaserPhotonEnergy,LaserNormIdist,LaserAbsCrossSection,LaserEmCrossSection,RefractiveIndex,GainMediumRefractiveIndex,Areas)
%Function used at each iteration to calculate the population in each energy level 
% Steady-state populations (from Edward Bernhardi's thesis)
% The input variables PumpPower and LaserPower are the sum of Fwd and Bwd
% powers

%Rates
R21 = 1/tau2;
R01 = (PumpPower.*PumpNormIdist.*PumpAbsCrossSection./PumpPhotonEnergy + LaserPower.*LaserNormIdist.*LaserAbsCrossSection./LaserPhotonEnergy)./Areas; %1E16 term is dividing by area...mesh element area???
R10 = (PumpPower.*PumpNormIdist.*PumpEmCrossSection./PumpPhotonEnergy + LaserPower.*LaserNormIdist.*LaserEmCrossSection./LaserPhotonEnergy)./Areas + 1/tau1;

%Auxiliar variables
A = W_ETU*(1+R01./R21);
B = R01 + R10;
C = -R01.*NT;

%Steady-state populations
N1 = (sqrt(B.^2-4*A.*C)-B)./(2*A);
N2 = (W_ETU*N1.^2)./R21;
N0 = NT - (N1+N2);

%Mesh points that are not TeO2 are not doped
for i=1:1:numel(RefractiveIndex(:,1))
    for j=1:1:numel(RefractiveIndex(1,:))
        if RefractiveIndex(i,j) ~= GainMediumRefractiveIndex
            N1(i,j)=0;
            N2(i,j)=0;
            N0(i,j)=0;
        end
    end
end

end