

function [yi] = extent01(sc, ej, ni0) %sc is equal to alpha in Tosun 

    ni = ni0 + sum(sc.*ej, 1);
    nt = sum(ni);
    yi = ni/nt;

end