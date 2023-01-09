addpath(genpath(cd))
run('InletConditions.m')

clc
clear
close all

load("Variables Storage\CO2HYD_ReactionProperties.mat")
%%
f1 = figure('WindowState', 'maximized', 'Color', 'w','Name',"Jacket Temperature & configuration");
ej0 = [1; 2; 4];

T = linspace(480, 700, 20);
P = linspace(40, 150, 20);
H2_CO2 = 2.5;

[Trange, Prange] = meshgrid(T, P); % T varies in columns and P in rows
XC1 = zeros(size(Trange));

for j = 1:size(Trange, 1)
    for i = 1:size(Trange, 1)       
        T = Trange(i,j);
        P = Prange(i,j);
        n0 = 100;     ni0  = [n0/(1+H2_CO2), 0, n0/(1+H2_CO2)*H2_CO2, 0, 0, 0];
        options = optimoptions('fsolve', 'algorithm','trust-region','FiniteDifferenceType', 'central', 'FunctionTolerance', 1e-10);   % supress display
        ej = fsolve(@(ej) eqCO2HYD06(ej, sc, T, P, ni0, Tc, Pc, Vc, wfac), ej0, options); %comp = 6; N20 = 0;
        ni = ni0 + sum(sc.*ej, 1);
        XC1(i, j) = (ni(4)) / sum(ni);
    end

end

subplot(2,2,1)
surf(Trange, Prange, XC1, EdgeColor="interp")
title('Sensitivity analysis on equilibrium function', Interpreter='latex')
xlabel('Temperature (K)', Interpreter='latex')
ylabel('Pressure (bar)', Interpreter='latex')
zlabel('Methanol molar fraction in equilibrium', Interpreter='latex')
cb = colorbar;                                     % create and label the colorbar
cb.Label.String = 'Methanol molar fraction';
view(45, 45)



%%
T = linspace(480, 700, 20);
H2_CO2 = linspace(0.5, 5, 20);
P = 150;

[Trange, H2_CO2range] = meshgrid(T, H2_CO2);
XC2 = zeros(size(Trange));

for j = 1:size(Trange, 1)
    for i = 1:size(Trange, 1)       
        T = Trange(i,j);
        H2_CO2 = H2_CO2range(i,j);
        n0 = 100;     ni0  = [n0/(1+H2_CO2), 0, n0/(1+H2_CO2)*H2_CO2, 0, 0, 0];
        options = optimoptions('fsolve', 'algorithm','trust-region','FiniteDifferenceType', 'central', 'FunctionTolerance', 1e-10);   % supress display
        ej = fsolve(@(ej) eqCO2HYD06(ej, sc, T, P, ni0, Tc, Pc, Vc, wfac), ej0, options); %comp = 6; N20 = 0;
        ni = ni0 + sum(sc.*ej, 1);
        XC2(i, j) = (ni(4)) / sum(ni);
    end
end

subplot(2,2,2)
surf(Trange, H2_CO2range, XC2, EdgeColor="interp")
title('Sensitivity analysis on equilibrium function', Interpreter='latex')
xlabel('Temperature (K)', Interpreter='latex')
ylabel('H2/CO2 ratio ', Interpreter='latex')
zlabel('Methanol molar fraction in equilibrium', Interpreter='latex')
cb = colorbar;                                     % create and label the colorbar
cb.Label.String = 'Methanol molar fraction';
view(45, 45)


%%

P = linspace(40, 140, 20);
H2_CO2 = linspace(0.5, 5, 20);
T = 480;

[Prange, H2_CO2range] = meshgrid(P, H2_CO2);
XC3 = zeros(size(Trange));

for j = 1:size(Trange, 1)
    for i = 1:size(Trange, 1)       
        P = Prange(i,j);
        H2_CO2 = H2_CO2range(i,j);
        n0 = 100;     ni0  = [n0/(1+H2_CO2), 0, n0/(1+H2_CO2)*H2_CO2, 0, 0, 0];
        options = optimoptions('fsolve', 'algorithm','trust-region','FiniteDifferenceType', 'central', 'FunctionTolerance', 1e-10);   % supress display
        ej = fsolve(@(ej) eqCO2HYD06(ej, sc, T, P, ni0, Tc, Pc, Vc, wfac), ej0, options); %comp = 6; N20 = 0;
        ni = ni0 + sum(sc.*ej, 1);
        XC3(i, j) =(ni(4)) / sum(ni);
    end
end

subplot(2,2,3)
surf(Prange, H2_CO2range, XC3, EdgeColor="interp")
title('Sensitivity analysis on equilibrium function', Interpreter='latex')
xlabel('Pressure (bar)', Interpreter='latex')
ylabel('H2/CO2 ratio ', Interpreter='latex')
zlabel('Methanol molar fraction in equilibrium', Interpreter='latex')
cb = colorbar;                                     % create and label the colorbar
cb.Label.String = 'Methanol molar fraction';

npath = fullfile(cd,'\Results\Surface Plot of Equilibrium.fig');
saveas(gcf, npath)
