addpath(genpath(cd))
run('InletConditions.m')

clc
clear
close all

load("Variables Storage\CO2HYD_ReactionProperties.mat")
%%
fsopt = optimoptions('fsolve', 'algorithm','trust-region','FiniteDifferenceType', 'central', 'FunctionTolerance', 1e-10, 'Display','off');   % supress display

n0 = 100;  %mol
ej0 = [1; 2;7];

A = []; %Variables matrix inequiality config     x1, x2, ... xn <= b
b = []; %Inequialities valor
Aeq = []; %Variables matrix equiality config
beq = []; %Equialities valor
ymin = [480, 50, 0.5]; %independent values min range search
ymax = [700, 150, 5];  %independent values max range search

y0 = [480, 150, 2.7];
Sol = fmincon(@(y) -eqCO2HYD07(y, n0, Tc, Pc, Vc, wfac, ej0, sc, fsopt), y0, A, b, Aeq, beq, ymin, ymax);

Topt = Sol(1)         %K
Popt = Sol(2)         %Bar
H2_CO2opt = Sol(3)    %-

XCmax = eqCO2HYD07(Sol, n0, Tc, Pc, Vc, wfac, ej0, sc, fsopt)
