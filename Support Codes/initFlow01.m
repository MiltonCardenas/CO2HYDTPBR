
%{
    R should have the same units as P, T, F
    Mi & Fi should have the same mole units
    Mi & Yi should be expressed as column vectors
    Mass units depends on Mi mass units
    Volumetric units depends on R volumetric units
    u units depends on cubic root of volumetric units and mole flow units
    Time units depends on mole Flow time units
%}

function [F, yi, V, u, rhog] = initFlow01(T, P, Fi, Ac, Mi, R)

    F  = sum(Fi);
    yi =  Fi/F;
    V   = F*T*R/P;
    rhog = sum(Mi.*yi) * P/(R*T);
    u = V/Ac;

end