
%formulation from fogler & Word Document
function [dPdW] = ergunW01(rhog, G, miug, rhocat, phi, At, dp, nt)

    dPdW = - G/(rhog*dp) * (1-phi)/(phi^3) * (150*(1-phi)*miug/dp + 1.75*G) * (1/((1-phi)*At*rhocat*nt));

    %kPa
end