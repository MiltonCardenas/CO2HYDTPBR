
%{

This equations assume equilibrium constant as
ln Kp° = 1/RT * (a1 + a2*T +a3*T^2 + a4*T^3 + a5*T^4 + a6*T^5 + a7*T*ln(T)

The input is a vector with dimensions (ix7) which contains the coefficients for the polynom below
where i is the number of reactions. The output is a vector which contains the calculated Ki variables 

R constant is implemented in the function and its value is R = 8.314; %J / mol K
Be careful about the R units and de T and coefficent units
%}

function [K1, K2, K3] = Keq_CO2HYD01(T)

R = 8.314472;
K = zeros(1, 3);

c = [7.44140e4, 1.89260e2, 3.2443e-2, 7.0432e-6, -5.6053e-9, 1.0344e-12, -6.4364e1;
    -3.94121e4, -5.41516e1, -5.5642e-2, 2.5760e-5, -7.6594e-9, 1.0161e-12, 1.8429e1];

Tpol = [1, T, T^2, T^3, T^4, T^5, T*log(T)];

for i = 1:2
    K(i) = exp((R*T)^-1*sum(c(i,:).*Tpol));
end

K(3) = K(2)*K(1);

K1 = K(1);
K2 = K(2);
K3 = K(3);

end