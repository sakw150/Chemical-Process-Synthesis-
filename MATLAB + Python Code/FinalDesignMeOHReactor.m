%% Appendix C: Methanol Reactor Modeling  
clear 
clc 
clear figure 
%% Methanol Reaction Kinetics 
% r1: CO+2H2=CH3OH
% r2: CO2+H2=CO+H2O
% r3: CO2+3H2=CH3OH+H2O

%% Volume Domain 
volume_domain=linspace(0,2000,2000);
%% Initial Conditions 
F_CO= 537809.040  %kmol/hr
F_H2= 1363146.62  %kmol/hr 
F_CH3OH= 241.87   %kmol/hr 
F_CO2= 602.74     %kmol/hr 
F_H2O= 535364.100 %kmol/hr
T_0=582.1500      %K
P_0=80            %Bar
T_C0=621.1500

IC=[F_CO, F_H2, F_CH3OH, F_CO2, F_H2O, T_0, P_0, T_C0];

%% Solve ODE
[Xsol, Ysol]=ode45('FinalDesignMeOHReactorODE', volume_domain, IC)
%% Data Handling
F_CO=Ysol(:,1);
F_H2=Ysol(:,2);
F_CH3OH=Ysol(:,3);
F_CO2=Ysol(:,4);
F_H2O=Ysol(:,5);
T=Ysol(:,6);
P=Ysol(:,7);
T_C=Ysol(:,8);

%% Figure Plotting
fig1=figure(1); 
set(fig1,'Name',"Flow Rate");
plot(Xsol,F_CO)
hold on 
plot(Xsol,F_H2)
hold on 
plot(Xsol,F_CH3OH)
hold on 
plot(Xsol,F_CO2)
hold on 
plot(Xsol,F_H2O)
hold off 
xlabel('Volume (L)'); 
ylabel('Flow Rate (mol/hr)'); 
legend("F_{CO}","F_{H2}","F_{CH3OH}","F_{CO2}","F_{H2O}")
title('Component Flow Rates')

fig2=figure(2);
set(fig2,'Name','Temperature');
plot(Xsol,T)
xlabel('Volume (L)');
ylabel('Temperature(K)');
title('Temperature Energy Balance')

fig3=figure(3);
set(fig3,'Name','Coolant Temperature');
plot(Xsol,T_C)
xlabel('Volume (L)');
ylabel('Temperature (K)');
title('Coolant Temperature')


fig4=figure(4)
set(fig4,'Name','Pressure')
plot(Xsol,P)
xlabel('Volume (L)');
ylabel('Pressure (bar)');
title('Pressure')
