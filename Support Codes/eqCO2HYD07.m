
function [XC] = eqCO2HYD07(y, n0, Tc, Pc, Vc, w, ej0, sc, fsopt)

 T = y(1);
 P = y(2);
 H2_CO2 = y(3);

 ni0  = [n0/(1+H2_CO2), 0, n0/(1+H2_CO2)*H2_CO2, 0, 0, 0];
 ej = fsolve(@(ej) eqCO2HYD06(ej, sc, T, P, ni0, Tc, Pc, Vc, w), ej0, fsopt); %comp = 6; N20 = 0;
 ni = ni0 + sum(sc.*ej, 1);
 XC = (ni(4) / (ni0(1)+ ni0(3)))/(1/4) ;

end