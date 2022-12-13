%% Appendix C: Zinc Oxide Reactor Modeling 
clear 
clc 
clear figure 

%% ZnO Reaction Kinetics 
% r1: ZnO + H2S = ZnS + H2O

%% Volume Domain 
volume_domain=linspace(0,1,200);

%% Initial Conditions 
F_H2S=0.311; %kmol/hr
F_ZnO=1200.8; %kmol/hr %we could assume PFR-like behavior
F_ZnS=0; %kmol/hr 
F_H2O=0.69778e+06; %kmol/hr 
T_0=298; %K
P_0=35; %Bar
%T_C0=600

IC=[F_H2S, F_ZnO, F_ZnS, F_H2O, T_0, P_0 ];

%% Solve ODE
[Xsol, Ysol]=ode45('FinalDesignZnOReactorODE', volume_domain, IC);

%% Data Handling
F_H2S=Ysol(:,1);
F_ZnO=Ysol(:,2);
F_ZnS=Ysol(:,3);
F_H2O=Ysol(:,4);
T=Ysol(:,5);
P=Ysol(:,6);
%% Figure Plotting

fig1=figure(1); 
set(fig1,'Name',"Concentration of H2S");
plot(Xsol,F_H2S)
hold on 
plot(Xsol,F_ZnS)
hold on  
xlabel('Volume (L)');  
ylabel('Flow Rate (kmol/hr)'); 
legend("F_{H2S}","F_{ZnS}")
title('Component Flow Rates')

fig2=figure(2);
set(fig2,'Name','Temperature');
plot(Xsol,T)
xlabel('Volume (L)');
ylabel('Temperature(K)');
title('Temperature Energy Balance')

fig4=figure(4);
set(fig4,'Name','Pressure')
plot(Xsol,P)
xlabel('Volume (L)');
ylabel('Pressure (bar)');
title('Pressure')
