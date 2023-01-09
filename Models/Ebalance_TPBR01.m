
function [dTdw] = Ebalance_TPBR01(T, TJ, F, Fi, r, U, As, hrxnj)
    
    r = r';
    cpmg = CpG_CO2HYD01(T, Fi);
    dTdw = (U*As*(TJ-T) - hrxnj(1)*r(1) - hrxnj(2)*r(2) - hrxnj(3)*r(3)) / (F*cpmg);
   
end