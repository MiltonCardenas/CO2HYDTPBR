
% X = Ej

function X = eqCO2HYD06(ej, sc, T, P, ni0, Tc, Pc, Vc, w)
    
    comp = 6;
    ni = ni0 + sum(sc.*ej, 1);
    nt = sum(ni);
    yi = ni/nt;
    
    [~, ~, fic] = ZSRK(comp, yi, Tc, Pc, Vc, w, T, P); 
    
    fic_CO2   = fic(1);
    fic_CO    = fic(2);
    fic_H2    = fic(3);
    fic_CH3OH = fic(4);
    fic_H2O   = fic(5);
    
    [Kp1, Kp2, Kp3] = Keq_CO2HYD01(T);
    X(1) = Kp1*(fic_CO)*((fic_H2)^2) - fic_CH3OH;
    X(2) = Kp2*(fic_CO2)*(fic_H2)-(fic_CO)*(fic_H2O);
    X(3) = Kp3*fic_CO2*(fic_H2^3)-fic_CH3OH*fic_H2O;

end