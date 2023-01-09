

%the K units depends on kjref units

function [k1, k2, k3, kCO, kCO2, kH2O_H2] = k_CO2HYD02(T)
    
    R = 8.314;

    k1 = 2.69e7*exp(-109900/(R*T));
    k2 = 7.31e8*exp(-123400/(R*T));
    k3 = 4.36e2*exp(-65200/(R*T));
    
    kCO = 7.99e-7*exp(58100/(R*T));
    kCO2 = 1.02e-7*exp(67400/(R*T));
    kH2O_H2  = 4.13e-11*exp(104500/(R*T));

    
end