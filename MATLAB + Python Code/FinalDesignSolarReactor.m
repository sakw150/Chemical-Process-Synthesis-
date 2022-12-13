%% Appendix C: Solar Reactor Modeling 
clear 
clc 
clear figure 
%% Solar Reaction Kinetics 
% r1: C6H1005+H2O=6CO+6H2
% r2: C10H1203+7H20=10CO+13H2
% r3: CO+H2O=CO2+H2 REVERSIBLE
% r4: CH4+H20=CO+3H2

%% Volume Domain 
volume_domain=linspace(0,50,300);
%% Initial Conditions 
F_C6H10O5= 13.414*0.2778  %kmol/hr
F_C10H1203= 37.7349*0.2778  %kmol/hr 
F_H2O= 266380*0.2778   %mol/s
F_H2= 5.57321*10^6*0.2778     %kmol/hr 
F_CO= 2.78661*10^6*0.2778 %kmol/hr
F_CO2 = 0*0.2778 %kmol/hr
F_CH4 = 259267*0.2778 %kmol/hr
T_0=298        %K


IC=[F_C6H10O5, F_C10H1203, F_H2O, F_H2, F_CO, F_CO2, F_CH4, T_0];

%% Solve ODE
[Xsol, Ysol]=ode45('FinalDesignSolarReactorODE', volume_domain, IC)

%% Data Handling
F_C6H10O5=Ysol(:,1);
F_C10H1203=Ysol(:,2);
F_H2O=Ysol(:,3);
F_H2=Ysol(:,4);
F_CO=Ysol(:,5);
F_CO2=Ysol(:,6);
F_CH4=Ysol(:,7);
T=Ysol(:,8);


%% Figure Plotting (Need To Change) 

fig1=figure(1); 
set(fig1,'Name',"Flow Rate");
plot(Xsol,F_C6H10O5)
hold on 
plot(Xsol,F_C10H1203)
hold on 
plot(Xsol,F_H2O)
hold on 
plot(Xsol,F_H2)
hold on 
plot(Xsol,F_CO)
hold on 
plot(Xsol, F_CO2)
hold on 
plot(Xsol, F_CH4)
hold off 
xlabel('Volume (L)'); 
ylabel('Flow Rate (mol/hr)'); 
legend("F_{C6H10O5}","F_{C10H1203}","F_{H2O}","F_{H2}","F_{CO}", "F_{CO2}", "F_{CH4}")
title('Component Flow Rates')

fig2=figure(2);
set(fig2,'Name','Temperature');
plot(Xsol,T)
xlabel('Volume (L)');
ylabel('Temperature(K)');
title('Temperature Energy Balance')

