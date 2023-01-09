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
%%   SINGLE CASE - VARIABLES PROFILES IN A TPBR REACTOR @(Inlet conditions, Reactor properties, Reaction properties)

Jconfig2 = "countercurrent";
[X, Y] = ode15s( @(w, y) ONEDS_CO2HYDTPBR01(w, y, T0, P0, F0, V0, rhog0, miug, U, FJ, CpJ, As, hrxnj, rhocat, phi, At, L, nt, dp, sc, Jconfig2), wspan, y0);
XC = Y(:, 4)/F0/(1/4);

ej0 = [1;2;4];
a = size(Y,1);
options = optimoptions('fsolve', 'algorithm','trust-region','FiniteDifferenceType', 'central', 'FunctionTolerance', 1e-10);   % supress display
ej = fsolve(@(ej) eqCO2HYD06(ej, sc, TJ0, P0/100, Fi0', Tc, Pc, Vc, wfac), ej0, options); %comp = 6; N20 = 0;
Fieq = Fi0' + sum(sc.*ej, 1);
Xeq = (Fieq(4)) / (Fi0(1) + Fi0(3))/(1/4);
Xeq = ones(1,a)*Xeq;

eq_reach = XC(end)/Xeq(end) *100;
%% GRAPHICS
f1 = figure('WindowState', 'maximized', 'Color', 'w','Name',"SINGLE CASE - VARIABLES PROFILES IN A TPBR REACTOR");
left_color =[0 0 0]; right_color =[0.4660 0.6740 0.1880]; set(f1,'defaultAxesColorOrder',[left_color; right_color])
lynestyle = [":", ":", ":", "-", ":", ":"];
colorstyle = ["#000000", '#00841a',"#000000",'#00841a','#00841a', "#000000"]; 
markers =  ["o", "v", "*", "x", "s", "^"]; 
legends = ["CO2", "CO", "H2", "CH3OH", "H2O", "N2"];

subplot(2,2,1)
n = length(Y(:, 1));

gs = [];
for i = 1:6
    if i == 4 || i == 5 || i == 2
        yyaxis right
        g = plot(X, Y(:,i), Color=colorstyle(i), Marker=markers(i), LineStyle=lynestyle(i), MarkerIndices=round(linspace(1,n,10),0), MarkerSize=5,LineWidth=0.7);
    else 
        yyaxis left
        g = plot(X, Y(:,i), Color=colorstyle(i), Marker=markers(i), LineStyle=lynestyle(i), MarkerIndices=round(linspace(1,n,10),0), MarkerSize=5,LineWidth=0.5);    
    end
    gs(i) = g;
    hold on
end

title('Species profile in the reactor', Interpreter="latex")
xlim([0, wcat])
xlabel('Catalyst (Kg)', Interpreter='latex')
ylabel('Molar flow (mol/s)', Interpreter='latex')
lgd1 = legend(gs, ["CO2", "CO", "H2", "CH3OH", "H2O", "N2"], Interpreter='Latex');
lgd1.Title.Interpreter = 'Latex';
yyaxis right
ylabel('Molar flow (mol/s)', Interpreter="latex")
hold off

subplot(2,2,2)
a = plot(X,Y(:,7), Color='k');
hold on
b = plot(X,Y(:,9), Color='r');
hold off

title('Temperature profile in the reactor', Interpreter="latex")
xlim([0, wcat])
xlabel('Catalyst (Kg)', Interpreter='latex')
ylabel('Temperature (K)', Interpreter='latex')
legend(["Reactor Temperature", "Jacket Temperature"], Interpreter='Latex')

subplot(2,2,3)
plot(X,Y(:,8)/100, 'Color', '#000000')
title('Pressure profile in the reactor', Interpreter="latex")
xlim([0, wcat])
xlabel('Catalyst (Kg)', Interpreter='latex')
ylabel('Pressure (bar)', Interpreter='latex')
hold off

subplot(2, 2, 4)
plot(X, XC, 'Color', '#000000')
hold on
plot(X, Xeq, '--b')
title(['Equilibrium reach is ', num2str(eq_reach), ' \%'], Interpreter="latex")
xlim([0, wcat])
xlabel('Catalyst (Kg)', Interpreter='latex')
ylabel('Conversion', Interpreter='latex')
lgd2 = legend("Reactor Conversion", "Equilibrium Conversion - ${{\rm{X}}_{\rm{C}}} = \frac{{{{\left( {{{\rm{n}}_{{\rm{C}}{{\rm{H}}_3}{\rm{OH}}}}} \right)}_{{\rm{eq}}}}/\left( {{{\left( {{{\rm{n}}_{{\rm{C}}{{\rm{O}}_2}}}} \right)}_0} + {{\left( {{{\rm{n}}_{{{\rm{H}}_2}}}} \right)}_0}} \right){\rm{\;}}}}{{{{\rm{\alpha }}_{{\rm{C}}{{\rm{H}}_3}{\rm{OH}}}}/\left( {{{\rm{\alpha }}_{{\rm{C}}{{\rm{O}}_2}}} + {{\rm{\alpha }}_{{{\rm{H}}_2}}}} \right){\rm{\;}}}}{\rm{\;}}\left( {{\rm{eq}}.16} \right)$",Interpreter='Latex');
lgd2.Location = 'best';
hold off

%-----------------------------
npath = fullfile(cd,'\Results\CO2HYD TPBR.fig');
saveas(gcf, npath)
