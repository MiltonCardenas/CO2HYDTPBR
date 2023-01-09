
function [FCH3OH] = kin_CO2HYD02(Y, sc)

T = Y(1);
P = Y(2);
yi = Y(3:end);

Pi = P*yi;

pCO2 =   Pi(1);
pCO =    Pi(2);
pH2 =    Pi(3);
pCH3OH = Pi(4);
pH2O =   Pi(5);
pN2 =    Pi(6);

% KCO, KGWS, KCO2
[K1 , K2, K3] = Keq_CO2HYD01(T);
[k1, k2, k3, kCO, kCO2, kH2O_H2] = k_CO2HYD02(T);            %  Ej y R are the same units?

% CO hydrog     (r1)
r1 = k1*kCO*(pCO*(pH2^1.5) - pCH3OH/((pH2^0.5)*K1)) / ((1+kCO*pCO+kCO2*pCO2)*((pH2^0.5)+kH2O_H2*pH2O));
% RGWS          (r2)
r2 = k2*kCO2*(pCO2*pH2 - pH2O*pCO/K2) / ((1+kCO*pCO+kCO2*pCO2)*((pH2^0.5)+kH2O_H2*pH2O));
% CO2 hydrog    (r2)
r3 = k3*kCO2*(pCO2*(pH2^1.5) - pCH3OH*pH2O/((pH2^1.5)*K3)) / ((1+kCO*pCO+kCO2*pCO2)*((pH2^0.5)+kH2O_H2*pH2O));

r = [r1, r2, r3]';

ri(1:6) = sum(r.*sc);
FCH3OH = ri(4);

end