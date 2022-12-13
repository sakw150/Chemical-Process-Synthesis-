function [per_ode]=ode(X,Y)
%% Output Variables 
F_C6H10O5=Y(1);  
F_C10H1203=Y(2);
F_H2O=Y(3);
F_H2=Y(4);
F_CO=Y(5);
F_CO2=Y(6);
F_CH4=Y(7);
T=Y(8);
%% Internal Variables 
R=8.314*0.001   %kJ/(mol*K)
X=0.90 % conversion
phi = 4.7 %geometric factor

P=35


K1=(2.51*10^3)*exp(-112.6/(R*T))
K2=(6.74*10^-2)*exp(-37.3/(R*T))
K3=(3.04*10^-1)*exp(-36.6/(R*T))

H_rxn1 = 336.791;   %kJ/mol cellulose
H_rxn2 = 2459.472;    %kJ/mol lignin
H_rxn3 = 175.181 ;   %kJ/mol methane
H_rxn4 = -210.905; %kJ/mol water-gas reaction


Cp1 =-0.011704*T+0.00067207*((T^2)/2) % Cellulose
Cp2 = 0.0314317*T+0.0003944*((T^2)/2) % Lignin
Cp3 = (0.033363)+0.02679*(0.00026105/T/sinh(0.00026105/T))^2+0.008896*(0.001169/T/cosh(0.001169/T))^2 % H20
Cp4 = (0.027617)+0.00956*(0.002466/T/sinh(0.002466/T))^2+0.00376*(0.0005676/T/cosh(0.0005676/T))^2% H2
Cp5 = (0.029108)+0.008773*(0.0030851/T/sinh(0.0030851/T))^2+0.008455*(0.0015382/T/cosh(0.0015382/T))^2 %CO
Cp6 = (0.02937)+0.03454*(0.001428/T/sinh(0.001428/T))^2+0.0264*(0.000588/T/cosh(0.000588/T))^2 % CO2
Cp7 = (0.033298)+0.079933*(0.0020869/T/sinh(0.0020869/T))^2+0.041602*(0.00099196/T/cosh(0.00099196/T))^2 % CH4



h_w = 0.0122 % W/m^2K 
u_g = 3.2*10^-5 % Pa*s
A_i = 0.01824 % m^2
V_m = 0.365 % m^3
p_c = 15.2701 %kg/m^3
d_i = 0.1524 %m
h_m = 0.7189 %W/m^2K
T_m = 1723.15 %K
MW_c = 0.17117255 %kg/mol
T_w = 1723.15 %K
A_m = 109.5 %m^3
F_m = 51.1489*(1000/3600) %kmol/hr

F_T0 = (8.64498*10^6)*(1000/3600); % Sum of all entering flow rates

F_tot=F_C6H10O5+F_C10H1203+F_H2O+F_H2+F_CO+F_CO2+F_CH4
P_C6H10O5 = F_C6H10O5/F_tot * P
P_C10H1203 = F_C10H1203/F_tot * P
P_H2O = F_H2O/F_tot * P
P_H2 = F_H2/F_tot * P
P_CO = F_CO/F_tot * P
P_CO2 = F_CO2/F_tot * P
P_CH4 = F_CH4/F_tot * P

k1=(K1*P_H2O)/(1+K2*P_H2O+K3*P_H2) %kJ/mol*bar*s
k2=(2.78*10^3)*exp((-1.26*10)/(R*T))
k3=(-9.59*10^4)*exp((-4.66*10)/(R*T))
k4=(3.0*10^8)*exp((-1.3*10^2)/(R*T))
%% ODE Definitions 
r1=k1*(1-X)*sqrt(1-phi*log(1-X))
r2=k2*P_CO*P_H2O
r3=k3*P_CO2*P_H2
r4=k4*P_CH4*P_H2O


dF_C6H10O5=-r1
dF_C10H1203=-r1
dF_H2O=-r1-r2+r3-r4
dF_CO=r1-r2+r3+r4
dF_H2=-r1+r2-r3+3*r4
dF_CO2=r2-r3
dF_CH4=-r4


Hpg=((h_m*A_m*F_m*MW_c)/(u_g*V_m*p_c))*(T-T_m)
Hwg=(h_w*pi*d_i)*(T-T_w)
SumFiCpi=(F_C6H10O5*Cp1)+(F_C10H1203*Cp2)+(F_H2O*Cp3)+(F_H2*Cp4)+(F_CO*Cp5)+(F_CO2*Cp6)+(F_CH4*Cp7)

dT= (A_i*((H_rxn1*r1)+(H_rxn2*r1)+(H_rxn3*r4)+(H_rxn4*r2)+(H_rxn4*r3))-A_i*Hpg-A_i*Hwg)/SumFiCpi


per_ode=[dF_C6H10O5; dF_C10H1203; dF_H2O; dF_H2; dF_CO; dF_CO2; dF_CH4; dT] ;
end 








