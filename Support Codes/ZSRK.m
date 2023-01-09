
function [Z, fugCoeff, fug] = ZSRK(comp, y, Tc, Pc, Vc, w, T, P) %in function of Z is not necessary R

Tr = T./Tc; %system temperature & pressure or species pressure
Pr = P./Pc;
alpha = (1+(0.480+1.574.*w-0.176.*(w.^2)).*(1-sqrt(Tr))).^2;

A = 0.42748*(Pr./(Tr.^2)).*alpha;
B = 0.08664*(Pr./Tr);

kij = zeros(comp);
Aij = zeros(comp);
amix = zeros(1,comp);

% estimation of binary interaction parameter
for i=1:comp
    num = Vc(i)^(1/6) * Vc.^(1/6);
    denom =  Vc(i)^(1/3) +  Vc.^(1/3);

    kij(i,:) = 1 - 8.*(num./denom).^3;
    kij(i,i) = 0;

    Aij(i,:) = (1-kij(i,:)).*sqrt(A(i).*A);
    amix(i) = sum(y(i).*y.*Aij(i,:));
end

A  = sum(amix);  
B = sum(y.*B);

p = -1; 
q = A-B-(B^2); 
r = -A*B;

Z0 = roots([1, p, q, r]);

Z = Z0(imag(Z0)==0);
Z = Z(1); % iff necessary, choose the higher Z root

fugCoeff = zeros(1,comp);
for i = 1:comp
    fugCoeff(i) = exp(B./B.*(Z-1) - log(Z-B) - A/B.*(2*sum(y.*Aij(i,:))/A - B./B).*log(1+(B./Z)));
end
fug = fugCoeff.*P.*y;
end