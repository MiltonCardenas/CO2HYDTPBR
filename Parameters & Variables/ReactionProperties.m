
clc
clear

%% CO2HYD REACTION SYSTEM PROPERTIES
%STOICHIOMETRIC COEFFICENTS -> CO2 CO H2 CH3OH H2O N2
sc =    [0,-1, -2, 1, 0, 0;               % rCO
        -1, 1, -1, 0, 1, 0;               % rRGWS
        -1, 0, -3, 1, 1, 0];              % rCH3OH

% SRK EoS INPUT PARAMETERS -> CO2 CO H2 CH3OH H2O N2
Mi = [44.009, 28.01, 2.016, 32.042, 18.01, 28.014]'/1000;       % Kg/mol (molecular weight)
Tc = [304.179, 134.18, 33.18, 512.67, 647.10, 126.20 ];         % K      (critical temperature)
Pc = [73.8, 37.078, 12.929, 80.558, 220.718, 33.98];            % bar    (critical pressure)
Vc = [94.28, 90.064, 65.03, 117.57, 55.989, 90.1 ];             % m3/mol (critical volume)
wfac  = [0.22551, 0.041035, -0.22014, 0.56197, 0.34417, 0.037]; % -      (accentric factor)

% REACTION ENTHALPY PER MOL OF CO2 -> r1, r2, r3
hrxnj = [-90.5, 41.4, -49.4 ]*1000;  % J/mol (reaction enthalpy)

%REACTION EXTENT FIRST GUESS
ej0 = [1; 2; 3];

% ----------------------------------------
save(fullfile(cd,'\Variables Storage\CO2HYD_ReactionProperties.mat'))

