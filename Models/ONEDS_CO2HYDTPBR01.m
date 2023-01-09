
function [dydw] = ONEDS_CO2HYDTPBR01(w, y, T0, P0, F0, V0, rhog0, miug, U, FJ, CpJ, As, hrxnj, rhocat, phi, At, L, nt, dp, sc, varargin)

Jconfig = varargin{1};

Fi = y(1:6);
T = y(7);
P = y(8);
TJ = y(9);

F = sum(Fi);
yi = Fi/F;
Pi = yi*P;  %rectify

[G, rhog, ~]= vflow_TPBR01(T, P, F, T0, P0, F0, V0, rhog0, phi, At, nt, L); %stationary deduction
r = kin_CO2HYD02(T, Pi);

dydw(1:6) = sum(r.*sc);
dydw(7) = Ebalance_TPBR01(T, TJ, F, Fi, r, U, As, hrxnj); %temperature modelling
dydw(8) = ergunW01(rhog, G, miug, rhocat, phi, At, dp, nt);
dydw(9) = Ebalance_JacketTPBR01(T, TJ, FJ, CpJ, U, As, phi, rhocat, Jconfig);

dydw = dydw';

end