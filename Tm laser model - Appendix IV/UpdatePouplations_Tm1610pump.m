function [N0,N1,N3] = UpdatePouplations_Tm1610pump(tau1,tau3,W_ETU,PumpPower,LaserPower,NT,PumpAbsCrossSection, PumpEmCrossSection, PumpPhotonEnergy,LaserPhotonEnergy,LaserAbsCrossSection,LaserEmCrossSection,BranchingRatio,W_RETU,ModeArea_pump,ModeArea_laser,OverlapTellurite_pump,OverlapTellurite_laser)
%Function used at each iteration to calculate the population in each energy level steady-state populations (see Chapter 2)
%To help with the fsolve convergence, I wrote  the equations in terms of the normalized populations (Ni/NT). Thanks to
%Arthur Mendez-Rosales for this suggestion.
% The input variables PumpPower and LaserPower are the sum of Fwd and Bwd powers

%Rates
R31 = (1-BranchingRatio)/tau3;
R01 = (PumpPower.*OverlapTellurite_pump.*PumpAbsCrossSection./PumpPhotonEnergy)/ModeArea_pump + (LaserPower.*OverlapTellurite_laser.*LaserAbsCrossSection./LaserPhotonEnergy)./ModeArea_laser; 
R10 = (PumpPower.*OverlapTellurite_pump.*PumpEmCrossSection./PumpPhotonEnergy)/ModeArea_pump + (LaserPower.*OverlapTellurite_laser.*LaserEmCrossSection./LaserPhotonEnergy)./ModeArea_laser + 1/tau1;
R32 = BranchingRatio/tau3;
            options = optimset('Display','off');
            F = @(x) [(R01.*x(1)-R10.*x(2)+R31.*x(3)+2.*(W_RETU.*x(1).*x(3).*NT-W_ETU.*x(2).*x(2).*NT));
                     -(1/tau3).*x(3)-(W_RETU.*x(1).*x(3).*NT-W_ETU.*x(2).*x(2).*NT);
                      (x(1)+x(2)+x(3)-1)];
            %Solve rate equations
            x=fsolve(F,[0.5 0.5 0.1],options);
            N0 = x(1)*NT;
            N1 = x(2)*NT;
            N3 = x(3)*NT;

        end
