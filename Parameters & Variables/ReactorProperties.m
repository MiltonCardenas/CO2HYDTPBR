
clc
clear

%% TUBULAR PACKED BED REACTOR PROPERTIES
%DIMENSIONAL PARAMETERS
L = 8;                            % m      (tube axial length)
Dt = 0.058;                       % m      (tube diameter)
nt = 150;                         % -      (tube number)
At = pi*((Dt/2)^2);               % m2     (individual tube cross-section area)

%PACKED BED PROPERTIES
dp = 0.006;                       % m        (particle catalyst diameter)
phi = 0.5;                        % -        (bed void fraction)  (empty volume / total bed volume)
rhocat = 1175;                    % kg/m3    (catalyst density)
As = 4*nt / (Dt*(1-phi)*rhocat);  % m2/kgcat (specific surface area of heat exchange)
wcat  = (1-phi)*At*rhocat*L*nt;   % KgCat

%ENERGY EXCHANGE SYSTEM PROPERTIES
U = 10;                           % J/m2sK            (global heat transfer coefficient)
Jconfig = "countercurrent";       % or "co-current"   (jacket configuration)  

% -------------------------------------
save(fullfile(cd,'\Variables Storage\TPBR_ReactorProperties.mat'))

