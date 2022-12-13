function [per_ode]=ode(X,Y)
%% Output Variables 
F_CO=Y(1);  
F_H2=Y(2);
F_CH3OH=Y(3);
F_CO2=Y(4);
F_H2O=Y(5);
T=Y(6);
P=Y(7);
T_C=Y(8)

%% Internal Variables 
R=8.314      
k1=(4.89*10^7)*exp(-113000/(R*T))
k2=(9.64*10^11)*exp(-152900/(R*T))
k3=(1.09*10^5)*exp(-87500/(R*T))

b_CO=(2.16*10^-5)*exp(46800/(R*T))
b_CO2=(7.05*10^-7)*exp(61700/(R*T))
b_H=(6.37*10^-9)*exp(84000/(R*T))      % Let b_H represent the ratio of b_H2O/(b_h2)^0.5

K1=10^(5139/T-12.621)
K2=10^(-2073/T-2.029)
K3=10^(3066/T-10.592)

H_rxn1 = 162.3387404;     %kJ/mol 
H_rxn2 = -192.8416736;    %kJ/mol
H_rxn3 = -206.7025002 ;   %kJ/mol


Cp1 = (29108)+8773*(3085.1/T/sinh(3085.1/T))^2+8455*(1538.2/T/cosh(1538.2/T))^2;
Cp2 = (27617)+9560*(2466/T/sinh(2466/T))^2+3760*(576.6/T/cosh(576.6/T))^2;
Cp3 = (39252)+87900*(37.6/T/sinh(37.6/T))^2+53000*(896.7/T/cosh(896.7/T))^2;
Cp4 = (29370)+34540*(1428/T/sinh(1428/T))^2+26400*(588/T/cosh(588/T))^2;
Cp5 = (33363)+26790*(2610.5/T/sinh(2610.5/T))^2+8896*(1169/T/cosh(1169/T))^2;

Cp10 = 5.11 ;    % Assuming that coolant specific heat is constant  

U_catalyst = 368.3/273.15;

density=17.59 
BP=0.41  % Bed Porosity 
rho_bulk = density*(1-BP) ;
phi = 0.50; % Void Frac 
rho_intrinsic = 440;

D_I=0.0254 
tubes=16
L=0.48 % meters 
a_s=D_I*L*tubes
a_c=pi*(D_I/2)^2*tubes
m_c = 10; % flowrate of coolant 

G =(7.63648*10^6)/(a_s*1000)*density  % Superficial mass velocity 
g_c =1  ; 

D_p = 1/8*D_I; % Heuristic
mu = 0.5; 

beta_0 =((G*(1-phi))/(rho_intrinsic*g_c*D_p*phi^3))*((150*(1-phi)*mu)/D_p+1.75*G)


F_T0 = 2.437*10^6; % Sum of all entering flow rates

F_tot=F_CO+F_CO2+F_H2+F_H2O+F_CH3OH;
P_CO = (F_CO/F_tot * P);
P_CO2 = F_CO2/F_tot * P;
P_H2 = F_H2/F_tot * P;
P_H2O = F_H2O/F_tot * P;
P_CH3OH = F_CH3OH/F_tot * P;

%% ODE Definitions 
r1=k1*b_CO*((P_CO*P_H2^(3/2)-P_CH3OH/(P_H2^0.5*K1))/((1+b_CO*P_CO+b_CO2*P_CO2)*(P_H2^0.5+b_H*P_H2O)))
r2=k2*b_CO2*((P_CO*P_H2-P_CO*P_H2O/K2)/((1+b_CO*P_CO+b_CO2*P_CO2)*(P_H2^(0.5)+b_H*P_H2O)))
r3=k3*b_CO2*((P_CO2*P_H2^(3/2)-((P_CH3OH*P_H2O)/(P_H2^(3/2)*K3)))/((1+b_CO*P_CO+b_CO2*P_CO2)*(P_H2^0.5+(b_H*P_H2O))))


dF_CO=-r1+r2
dF_H2=-2*r1-r2-3*r3
dF_CH3OH=r1+r3
dF_CO2=-r2-r3
dF_H2O=r2+r3

dT=(((U_catalyst/rho_bulk)*(T_C-T)+r1*H_rxn1+r2*H_rxn2+r3*H_rxn3)/(F_CO*Cp1+F_H2*Cp2+F_CH3OH*Cp3+F_CO2*Cp4+F_H2O*Cp5))

dT_c=U_catalyst*a_s*(T-T_C)/(m_c*Cp10);

dP=-((beta_0/(a_c*(1-phi)))*((80)/P)*(T/T_C)*(F_tot/F_T0))

per_ode=[dF_CO; dF_H2; dF_CH3OH; dF_CO2; dF_H2O; dT; dP; dT_c] ;
end 


