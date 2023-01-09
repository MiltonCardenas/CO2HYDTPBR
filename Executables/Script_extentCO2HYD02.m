addpath(genpath(cd))
run('InletConditions.m')

clc
clear
close all

load("Variables Storage\CO2HYD_ReactionProperties.mat")
%%
n = 25;
T = linspace(480, 700, n);
P = linspace(40, 140, n);
H2_CO2 = linspace(1, 4, n);
ej0 = [1; 2; 7];

[Trange, Prange, H2_CO2range] = meshgrid(T,P,H2_CO2);
XC = zeros(size(Trange));

for j = 1:length(Trange)
    for i = 1:length(Trange)
        for s = 1:length(Trange)

            T = Trange(i, j, s);
            P = Prange(i, j, s);
            H2_CO2 = H2_CO2range(i,j,s);
            n0 = 100;     ni0  = [n0/(1+H2_CO2), 0, n0/(1+H2_CO2)*H2_CO2, 0, 0, 0];
            options = optimoptions('fsolve', 'algorithm','trust-region','FiniteDifferenceType', 'central', 'FunctionTolerance', 1e-10);   % supress display
            ej = fsolve(@(ej) eqCO2HYD06(ej, sc, T, P, ni0, Tc, Pc, Vc, wfac), ej0, options); %comp = 6; N20 = 0;
            ni = ni0 + sum(sc.*ej, 1);
            XC(i,j,s) = ni(4) / (ni0(1)+ni0(3));

        end
    end
end


points = 625;
TT = zeros(1, points);
PP = zeros(1, points);
HH = zeros(1, points);
XC1 = zeros(1, points);

for coun = 1:points
    
    a = 1+abs(round((n-1)*rand(1)));
    b = 1+abs(round((n-1)*rand(1)));
    c = 1+abs(round((n-1)*rand(1)));

    TT(coun) = Trange(a, b, c);
    PP(coun) = Prange(a, b, c);
    HH(coun) = H2_CO2range(a, b, c);

    XC1(coun) = XC(a,b,c) ;

end


f1 = figure('WindowState', 'maximized', 'Color', 'w','Name',"Jacket Temperature & configuration");
scatter3(TT,PP,HH,40, XC1,'filled')
ax = gca;
ax.XDir = 'reverse';
cb = colorbar;                                     % create and label the colorbar
cb.Label.String = 'Conversion';
title('Sensitivity analysis on equilibrium function', Interpreter='latex')
xlabel('Temperature (K)', Interpreter='latex')
ylabel('Pressure (bar)', Interpreter='latex')
zlabel('H2/CO2 ratio', Interpreter='latex')
view(45, 45)

npath = fullfile(cd,'\Results\Scatter Plot of Equilibrium.fig');
saveas(gcf, npath)

