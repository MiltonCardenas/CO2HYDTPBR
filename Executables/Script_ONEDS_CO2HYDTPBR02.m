
addpath(genpath(cd))
run('InletConditions.m')

clc
clear
close all
files = dir(fullfile(cd, 'Variables Storage'));           % Search the carpet in the current directory (cd) and get existing filenames and properties
for i = 1:length(files)
    if isfolder(files(i).name)
        nf = 1;
    else 
        load(files(i).name, '-mat')                               %only can read non-folder files
    end
end

%Fi0 is calculated in this code
%% Sensitivity analysis dependent variables range

H2CO20 = 1;             H2CO2f = 5;
P0 = 5000;              Pf = 15000;
T0 = 150 + 273.15;      Tf = 270 + 273.15;
n = 50;

%
H2_CO2range = linspace(H2CO20, H2CO2f, n);
Prange = linspace(P0, Pf, n);
Trange = linspace(T0, Tf, n);

Jconfig1 = "countercurrent";
Jconfig2 = "co-current";

ys_1 = zeros(n, 9);
ys_2 = zeros(n, 9);
ys_3 = zeros(n, 9);
ys_4 = zeros(n, 9);
ys_5 = zeros(n, 9);

XC1 = zeros(1, n);
XC2 = zeros(1, n);
XC3 = zeros(1, n);
XC4 = zeros(1, n);
XC5 = zeros(1, n);

i = 1;
%%
for H2_CO2 = H2_CO2range

    Fi0  = [F0/(1+H2_CO2), 0, F0/(1+H2_CO2)*H2_CO2, 0, 0, 0]';
    [F0, yi0, V0, u0, rhog0] = initFlow01(T0, Pf, Fi0, At, Mi, R);
    y0 = [Fi0; T0; Pf; TJ0];

    [X, Y] = ode15s(@(w, y) ONEDS_CO2HYDTPBR01(w, y, T0, P0, F0, V0, rhog0, miug, U, FJ, CpJ, As, hrxnj, rhocat, phi, At, L, nt, dp, sc, Jconfig), wspan, y0);
    XC1(i) = Y(end, 4)/(Fi0(1)+Fi0(3))/(1/4);
    ys_1(i, :) = Y(end, :);
    i = i+1;    
end

f1 = figure('WindowState', 'maximized', 'Color', 'w');
lynestyle = [":", ":", ":", "-", ":", ":"];
colorstyle = ["#000000", '#00841a',"#000000",'#00841a','#00841a', "#000000"]; 
markers =  ["o", "v", "*", "x", "s", "^"]; 
legends = ["CO2", "CO", "H2", "CH3OH", "H2O", "N2"];


subplot(2,2,1)
for i = 1:6
plot(H2_CO2range, ys_1(:, i), 'LineStyle', ':', 'Color', 'k', 'Marker', markers(i), 'MarkerIndices', round(linspace(1,n,10),0),'MarkerSize', 4);
hold on
end
legend(legends, Interpreter='latex')
title('Inlet H2/CO2 ratio sensitivity analysis', Interpreter='latex')
xlabel('Inlet H2/CO2 ratio', Interpreter='latex')
ylabel('Final molar flow (mol/s)', Interpreter='latex')

H2_CO2 = 2.6;
Fi0  = [F0/(1+H2_CO2), 0, F0/(1+H2_CO2)*H2_CO2, 0, 0, 0]';

i = 1;
for P0 = Prange
    
    [F0, yi0, V0, u0, rhog0] = initFlow01(T0, P0, Fi0, At, Mi, R);
    y0 = [Fi0; T0; P0; TJ0];
   
    [X, Y] = ode15s(@(w, y) ONEDS_CO2HYDTPBR01(w, y, T0, P0, F0, V0, rhog0, miug, U, FJ, CpJ, As, hrxnj, rhocat, phi, At, L, nt, dp, sc, Jconfig), wspan, y0);
    XC2(i) = Y(end, 4)/(Fi0(1)+Fi0(3))/(1/4);
    ys_2(i, :) = Y(end, :);
    i = i+1;
end

subplot(2,2,2)
for i = 1:6
plot(Prange/100, ys_2(:, i),'LineStyle', ':','Color', 'k', 'Marker', markers(i), 'MarkerIndices', round(linspace(1,n,10),0),'MarkerSize', 4);
hold on
end
title('Inlet pressure sensitivity analysis', Interpreter='latex')
xlabel('Inlet pressure (bar)', Interpreter='latex')
ylabel('Final molar flow (mol/s)', Interpreter='latex')
legend(legends, Interpreter='latex')

P0 = 5500;
i = 1;
for T0 = Trange
    
    [F0, yi0, V0, u0, rhog0] = initFlow01(T0, P0, Fi0, At, Mi, R);
    y0 = [Fi0; T0; P0; TJ0];
    %[X, Y] = ode15s(@(w, y) ONEDS_CO2HYDTPBR01(w, y, T0, P0, F0, V0, rhog0, miug, U, FJ, CpJ, As, hrxnj, rhocat, phi, At, L, nt, dp, sc, Jconfig1, Tc, Pc, Vc, wfac,comp, Zmix0), wspan, y0);
    [X, Y] = ode15s(@(w, y) ONEDS_CO2HYDTPBR01(w, y, T0, P0, F0, V0, rhog0, miug, U, FJ, CpJ, As, hrxnj, rhocat, phi, At, L, nt, dp, sc, Jconfig), wspan, y0);
    XC3(i) = Y(end, 4)/(Fi0(1)+Fi0(3))/(1/4);
    ys_3(i, :) = Y(end, :);
    i = i+1;
end

subplot(2,2,3)
for i = 1:6
plot(Trange, ys_3(:, i), 'LineStyle', ':', 'Color', 'k', 'Marker', markers(i), 'MarkerIndices', round(linspace(1,n,10),0),'MarkerSize', 4);
hold on
end
title('Inlet temperature sensitivity analysis', Interpreter='latex')
xlabel('Inlet temperature (K)', Interpreter='latex')
ylabel('Final molar flow (mol/s)', Interpreter='latex')
xlim([Trange(1), Trange(end)])
legend(legends, Interpreter='latex')

T0 = 150 + 273.15;
i=1;
for TJ0 = Trange

    [F0, yi0, V0, u0, rhog0] = initFlow01(T0, Pf, Fi0, At, Mi, R);
    y0 = [Fi0; T0; Pf; TJ0];

    % CounterCurrent
    [X, Y] = ode15s(@(w, y) ONEDS_CO2HYDTPBR01(w, y, T0, Pf, F0, V0, rhog0, miug, U, FJ, CpJ, As, hrxnj, rhocat, phi, At, L, nt, dp, sc, Jconfig), wspan, y0);
    XC4(i) = Y(end, 4)/(Fi0(1)+Fi0(3))/(1/4);

    % Co-Current
    [X1, Y1] = ode15s(@(w, y) ONEDS_CO2HYDTPBR01(w, y, T0, Pf, F0, V0, rhog0, miug, U, FJ, CpJ, As, hrxnj, rhocat, phi, At, L, nt, dp, sc, Jconfig2), wspan, y0);
    XC5(i) = Y1(end, 4)/(Fi0(1)+Fi0(3))/(1/4);

    ys_4(i, :) = Y(end, :); %Jconfig
    ys_5(i, :) = Y1(end, :);%Jconfig2
    i = i+1;
end

gs1 = [];
gs2 = [];

subplot(2,2,4)
for i = 1:6
 plot(Trange, ys_4(:, i), 'Color', 'r', 'LineStyle',':', 'Marker', markers(i), 'MarkerIndices', round(linspace(1,n,10),0),'MarkerSize', 4);

hold on
end

for i = 1:6
    plot(Trange, ys_5(:, i), 'Color', '#000000', 'LineStyle', '--', 'Marker', markers(i), 'MarkerIndices', round(linspace(1,n,10),0), 'MarkerSize', 4);
    hold on 
end

title('Jacket Temperature Sensitivity Analysis: Counter current(black), Co-Current(red)', Interpreter='latex')
xlabel('System Temperature (K)', Interpreter='latex')
ylabel('Final molar flow (mol/s)', Interpreter='latex')
xlim([Trange(1), Trange(end)])

npath = fullfile(cd,'\Results\Sensitivity Analysis in TPBR 01.fig');
saveas(gcf, npath)

%%

f2 = figure('WindowState', 'maximized', 'Color', 'w','Name',"Jacket Temperature & configuration");

subplot(2,2,1)
plot(H2_CO2range, XC1,'k')
title('Inlet H2/CO2 ratio effect on XC in TPBR', Interpreter='latex')
xlabel('H2/CO2 range ', Interpreter='latex')
ylabel('Feedstock conversion to methanol XC', Interpreter='latex')

subplot(2,2,2)
plot(Prange/100, XC2, 'k')
title('Inlet pressure effect on XC in TPBR', Interpreter='latex')
xlabel('Inlet pressure (bar)', Interpreter='latex')
ylabel('Feedstock conversion to methanol XC', Interpreter='latex')

subplot(2,2,3)
plot(Trange, XC3,'k--')
hold on
plot(Trange, XC4,'k')
hold on
plot(Trange, XC5,'r')
title('Inlet Temperature effect on XC in TPBR', Interpreter='latex')
xlabel('Inlet temperature (K)', Interpreter='latex')
ylabel('Feedstock conversion to methanol XC', Interpreter='latex')
xlim([Trange(1), Trange(end)])
lgd2 = legend('Inlet reactor temperature effect on XC', 'Inlet jacket temperature in countercurrent effect in XC','Inlet jacket temperature in co-current effect in XC', Interpreter='latex');
lgd2.Location= 'northwest';

npath = fullfile(cd,'\Results\Sensitivity Analysis in TPBR 02.fig');
saveas(gcf, npath)
