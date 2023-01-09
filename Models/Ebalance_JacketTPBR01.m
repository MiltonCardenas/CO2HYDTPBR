
function dTJdw = Ebalance_JacketTPBR01(T, TJ, FJ, CpJ, U, As, phi, rhocat, varargin)

Jconfig = varargin{1};
rhob = (1-phi)*rhocat;
  
if strcmp(Jconfig, "co-current") %compare strings function
    dTJdw = U*As/rhob * (T-TJ)/(FJ*CpJ);
elseif strcmp(Jconfig, "countercurrent")
    dTJdw = U*As/rhob * (TJ-T)/(FJ*CpJ);
else
    disp("Check spelling mistakes: only {co-current} or {countercurrent} available")
end

