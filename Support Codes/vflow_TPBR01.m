
function [G, rhog, V] = vflow_TPBR01(T, P, F, T0, P0, F0, V0, rhog0, phi, At, nt, L)


    V = V0 * F/F0 * T/T0 * P0/P;
    rhog = rhog0*V0/V;

    Ac = At*L*nt*phi;
    u = V/Ac;

    G = rhog*u;
end