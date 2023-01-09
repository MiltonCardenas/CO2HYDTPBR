
addpath(genpath(cd))
run('ReactorProperties.m'); run('ReactionProperties.m')

clc
clear
load('TPBR_ReactorProperties.mat', 'At', 'wcat'); load('CO2HYD_ReactionProperties.mat', 'Mi')

%% SPECIES INLET CONDITIONS
% CO2 CO H2 CH3OH H2O N2
T0 = 475;        % K     (Species inlet temperature)                 
P0 = 15000;      % kpa   (Species inlet pressure)
F0 = 100; 
H2_CO2 = 2.6;    % -     (Hydrogen Carbon dioxide inlet ratio) H2CO2 opt = estequiometric quantities,
Fi0  = [F0/(1+H2_CO2), 0, F0/(1+H2_CO2)*H2_CO2, 0, 0, 0]';   %mol/s  (species inlet molar flow)

R = 0.008314472;                                                  % kpam3/molK   kJ/molK  Gas constant  (rhoGas)
[F0, yi0, V0, u0, rhog0] = InitFlow(T0, P0, Fi0, At, Mi, R);      % mol/s, -, m3/s, m/s, kg/m3 -> initial total mole flow, initial mole fraction, initial volumetric flow, initial gas velocity, initial gas density 
miug = 2.584*10^-5;                                               % Pa*s  Gas viscocity

%Heat exchange system
TJ0 = 526;                              % K            Jacket Fluid temperature
FJ = 0.25    ;                          % mol/s
CpJ = 75.327;                           % cte T = 25°C

% Inlet Ode Conditions
y0 = [Fi0; T0; P0; TJ0];
wspan = [0, wcat];

% --------------------------
save(fullfile(cd,'\Variables Storage\InletConditions.mat'))

