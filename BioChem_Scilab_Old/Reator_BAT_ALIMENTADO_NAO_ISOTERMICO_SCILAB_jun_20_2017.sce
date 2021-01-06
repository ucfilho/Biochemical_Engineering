clc
clear
//codigo scilab sistema Edo
//Eng.Bioq prof.Ubirajara-fev 07 2014
function [f]=BACH(t,Var)
T1=Var(1);Tc=Var(2);Tj=Var(3);
V1=Var(4);X1=Var(5);S1=Var(6);P1=Var(7);
Pmax=Alfa*exp(Beta*T1);
MI_0=k0*exp(-E/(R*(T1+273)));
MI=MI_0*S1/(S1+Ks)*(1-P1/Pmax)^n*(1-X1/Xmax)^m;
if(V1 < Vdorna) F0=Falim; else F0=0;V1=Vdorna;end
LMDT=((T1-Tj)-(Tc-Tje))/log((T1-Tj)/(Tc-Tje)) ;
dT1_dt=(F0*T0-F1*T1-F0*T1+Fc*(Tc-T1)-deltaH*X1/(Rho*Cp*Yxs))/V1;
dTc_dt=Fc/Vc*(T1-Tc)-U*A/(Rho*Cp*Vc)*LMDT;
dTj_dt=Fj/Vj*(Tje-Tj)+U*A/(Rhoj*Cpj*Vj)*LMDT;
dV1_dt=F0;
dS1_dt=(F0*S0-F1*S1-F0*S1-V1*X1*MI/Yxs)/V1;
dX1_dt=(F0*X0-F1*X1-F0*X1+V1*X1*MI)/V1;
dP1_dt=(F0*P0-F1*P1-F0*P1+V1*X1*MI*(Yps/Yxs))/V1;
f(1)=dT1_dt;f(2)=dTc_dt;f(3)=dTj_dt;
f(4)=dV1_dt;f(5)=dX1_dt;f(6)=dS1_dt; f(7)=dP1_dt;
endfunction


// Programa principal
k0=4.50e10;E=1.54e4;R=1.987;
Yxs=0.033;Yps=0.445;Ks=1.6;Xmax=100;deltaH=-158;
n=3;m=0.9;a=-0.0676;A=4.5e10;
//Falim=100;
F1=0;Falim=80;T1=32;
Fc=180;Fj=180;
Rhoj=1000;Rho=950;
Cp=1;Cpj=1;
V0=50;Vc=20;Vj=20;Vdorna=300;
A=210;U=3500;
Alfa=895.3;Beta=-0.0676;
T1=32;Tc=30;Tj=30;T0=30;Tje=25;
V1=50;
X0=50;S0=180;P0=1;
t0=0;U0=[T1;Tc;Tj;V0;X0;S0;P0];
t=[1;1.5;2;2.5;3;3.5;4;4.5;5;5.5;6;6.5;7;7.5;8;8.5;9;9.5;10];
Solucao=ode(U0,t0,t,BACH);
disp("Time---------T1--------Tc-------------Tj-----------Vol--------Cel--------S----------P");
disp([t,Solucao']);
