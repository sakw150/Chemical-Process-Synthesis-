function [per_ode]=FinalDesignZnOReactorODE(X,Y)
%% Output Variables 
F_H2S=Y(1);  
F_ZnO=Y(2);
F_ZnS=Y(3);
F_H2O=Y(4);
T=Y(5);
P=Y(6);

%% Internal Variables 
R_0=0.0765; %um
Z_v=0.6086; %unitless
x=0.90; %unitless
D_eff = 9.2*10^-14; %m^2/s
k_1 = 0 ;%unitless
C_pg = 5.97635*10^-8; %kJ/hr
m = 1.29*10^8; %kg/hr
h_g = 420.90 ; %W/m^2 K
a_v = 1; %m
T_s = 483; %K
%unit conversions were added to each differential eqn

%For Ergun Equation
G = 12.14*10^5 ; %mol/cm^2 s
rho = 0.02489377 ; %mol/L
D_p = 0.153 ; %um
phi = 0.15 ; %unitless
mu = 0.000143 ; %kg m/s

% With changes in flow
F_tot=F_H2S+F_ZnO+F_ZnS+F_H2O;
P_H2S = F_H2S/F_tot * P;
P_ZnO = F_ZnO/F_tot * P;
P_ZnS = F_ZnS/F_tot * P;
P_H2O = F_H2O/F_tot * P;
R = 0.08314; %L bar/mol K   this is the gas constant used right below
C_H2S= P_H2S/R/T;

%% ODE Definitions 
r1= 4*pi*R_0^2*C_H2S*(R_0*((Z_v+(1-Z_v)*(1-x))^(1/3)-(1-x)^(1/3))/((Z_v+(1-Z_v)*(1-x))^(1/3)*(1-x)^(1/3)*D_eff)+1/(1-x)^(2/3)*k_1)^(-1)*10^20;


% Differentials Set Below 

dF_ZnO=-r1;
dF_H2S=-r1;
dF_ZnS=r1;
dF_H2O=r1;

dT=-h_g*a_v*(T-T_s)/(m*C_pg)*60^2*10^-3*34.1;

dP=(-G/(rho*D_p)*10^7)*((1-phi)/phi^3)*((150*(1-phi)*mu)/D_p*10^5/34.1+1.75*G)*10^-25; %some values are unit conversions

per_ode=[dF_ZnO; dF_H2S; dF_ZnS; dF_H2O; dT; dP] ;
end 


