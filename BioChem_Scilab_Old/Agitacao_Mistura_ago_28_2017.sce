clc
clear
Ro=1000;MI=1e-3;// m3/s kg/m.s
Vol=300;// m3
Np=0.28;//hidrofoil 3 pas
D=(4*Vol/%pi)^(1/3);disp(D,"D=");
Di=D/3;disp(Di,"Di=");
X=[300 30];//rpm


fs=1.25;
for j=1:2
    Nrpm=X(j);
    disp("----------------------------")
    disp(Nrpm,"Nrpm")
    Nrps=Nrpm/60;
    Re=Ro*Di^2*Nrps/MI;disp(Re,"Re=");
    Pot=Np*Ro*Di^5*Nrps^3;disp(Pot,"Pot=");
    i=1750/Nrpm;disp(i,"i=");
    Pmotor=fs*Pot/745.7;disp(Pmotor,"Pmotor");
    Torque=Pot/(2*%pi*Nrps);disp(Torque,"Torque");
    Tmotor=Torque/(i*fs);disp(Tmotor,"Tmotor");
end
