addpath(genpath(cd))
run('InletConditions.m')

clc
clear
close all

load("Variables Storage\CO2HYD_ReactionProperties.mat")

yi = [9.97504060381826e-08	0.0987814913096220	0.901217503593899	3.86715182727320e-07	3.03312749378339e-08	4.88299617155514e-07];
ni = 100*yi;
%ni = [0, 10, 30, 0, 0, 0];
%yi = ni/sum(ni);

n = 10;
T = linspace(480, 520, n);
P = linspace(30, 150, n);

[Trange, Prange] = meshgrid(T,P);

rCH3OH = zeros(size(Trange));
r1 = zeros(size(Trange));
r2 = zeros(size(Trange));
r3 = zeros(size(Trange));

for j = 1: size(Trange, 1)
    for i = 1:size(Trange, 1)
        
        T = Trange(i, j);
        P = Prange(i, j);

        [~,~,Pi] = ZSRK(6, yi, Tc, Pc, Vc, wfac, T, P);
        r = kin_CO2HYD02(T, Pi);

        r1(i,j) = r(1);
        r2(i,j) = r(2);
        r3(i,j) = r(3);
        rCH3OH(i, j) = sum(r.*sc(:, 4));

    end
end

f1 = figure('WindowState', 'maximized', 'Color', 'w','Name',"Sensitivity analysis of Kinetic function");
subplot(2,2,1)
surf(Trange, Prange, r1, EdgeColor="interp")
title('Sensitivity analysis on Kinetics function', Interpreter='latex')
xlabel('Temperature (K)', Interpreter='latex')
ylabel('Pressure ', Interpreter='latex')
zlabel('Reaction rate (r1) (mol/kgCat*s)', Interpreter='latex')
hold off


subplot(2,2,2)
surf(Trange, Prange, r2, EdgeColor="interp")
title('Sensitivity analysis on Kinetics function', Interpreter='latex')
xlabel('Temperature (K)', Interpreter='latex')
ylabel('Pressure ', Interpreter='latex')
zlabel('Reaction rate (r2) (mol/kgCat*s)', Interpreter='latex')

subplot(2,2,3)
surf(Trange, Prange, r3, EdgeColor="interp")
title('Sensitivity analysis on Kinetics function', Interpreter='latex')
xlabel('Temperature (K)', Interpreter='latex')
ylabel('Pressure ', Interpreter='latex')
zlabel('Reaction rate (r3) (mol/kgCat*s)', Interpreter='latex')

subplot(2,2,4)
surf(Trange, Prange, rCH3OH, EdgeColor="interp")
title('Sensitivity analysis on Kinetics function', Interpreter='latex')
xlabel('Temperature (K)', Interpreter='latex')
ylabel('Pressure ', Interpreter='latex')
zlabel('Reaction rate (rCH3OH) (mol/kgCat*s)', Interpreter='latex')
hold off

%%

T = linspace(480, 560, n);
P = linspace(30, 150, n);

[Trange, Prange] = meshgrid(T,P);
a = size(Trange);
r2 = zeros(size(Trange));
r3 = zeros(size(Trange));

for j = 1: size(Trange, 1)
    for i = 1:size(Trange, 1)
        
        T = Trange(i, j);
        P = Prange(i, j);

        [~,~,Pi] = ZSRK(6, yi, Tc, Pc, Vc, wfac, T, P);
        r = kin_CO2HYD02(T, Pi);

        r1(i, j) = r(1);
        r2(i, j) = r(2);
        r3(i, j) = r(3);
    end
end

npath = fullfile(cd,'\Results\Sensitivity analysis of Kinetic Function.fig');
saveas(gcf, npath)