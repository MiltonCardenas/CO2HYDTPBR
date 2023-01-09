addpath(genpath(cd))
run('InletConditions.m')

clc
clear
close all

load("Variables Storage\CO2HYD_ReactionProperties.mat")
%%   OPTIMIZATION CODE FOR CH3OH KINETCIS

options = optimoptions("fmincon","FiniteDifferenceType","central","MaxFunctionEvaluations",1e5,"MaxIterations",1e5);

A = []; %Variables matrix inequiality config     x1, x2, ... xn <= b
b = []; %Inequialities valor
Aeq = [0, 0, 1, 1, 1, 1, 1, 1; %sum yi = 1;
       0, 0, 0, 0, 0, 0, 0 , 0]; %Variables matrix equality  %limit reagent
beq = [1,0]; %Equialities valor

ymin = [480, 30, 0, 0, 0, 0, 0, 0]; %independent values min range search
ymax = [520, 150, 1, 1, 1, 1, 1, 1];  %independent values max range search
nonlin = []; %nonlinear cons

yi0 = [25, 0, 75, 0, 0, 0]/100;

y0 = [500, 100, yi0];
Sol = fmincon(@(Y) -kin_CO2HYD03(Y, sc), y0, A, b, Aeq, beq, ymin, ymax,nonlin, options);  %must return escalar value %IDEAL

Topt = Sol(1)
Popt = Sol(2)
yiOpt = Sol(3:8)

rCH3OH = kin_CO2HYD03(Sol, sc); %mols/kgCat*s
