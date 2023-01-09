
%{

This equations assume equilibrium constant as
ln Kp° = 1/RT * (a1 + a2*T +a3*T^2 + a4*T^3 + a5*T^4 + a6*T^5 + a7*T*ln(T)

The input is a vector with dimensions (ix7) which contains the coefficients for the polynom below
where i is the number of reactions. The output is a vector which contains the calculated Ki variables 

R constant is implemented in the function and its value is R = 8.314; %J / mol K
Be careful about the R units and de T and coefficent units
%}

function [K1, K2, K3] = Keq_CO2HYD02(T)

K2 = 10^(-2073/T + 2.029);
K3 = 10^(3066/T -10.592);

K1 = K3/K2;


end