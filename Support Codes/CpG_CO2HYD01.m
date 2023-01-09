
% Ideal gas formulation

function [cpg] =  CpG_CO2HYD01(T, Fi)

cpi_coeff = [29.268, -0.224, 2.653, -4.153, 20.057;
           29.651, -0.007, 0.183, -0.094, 0.108;
           27.004,  0.119, -0.241, 0.215, -0.615;
           36.155, -0.511, 2.215, -1.824, 4.899;
           33.763, -0.006, 0.224, -0.100, 0.110;
           29.802, -0.070, 0.174, -0.085, 0.093].* [1, 1e-1, 1e-4, 1e-7, 1e-11]; %ideal from tosun J/molK 200<T<800

    yi = Fi/sum(Fi);
    cpig = sum(cpi_coeff .* [1, T, T^2, T^3, T^4],2);  %sum rows
    cpg = sum(yi.*cpig);

end

