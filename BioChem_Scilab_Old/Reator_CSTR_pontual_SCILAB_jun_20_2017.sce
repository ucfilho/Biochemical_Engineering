
function dy=f(x,y)
    S=y(1);X=y(2);P=y(3);
    M=MI*S/(S+Ks)
    dy(1)=D*S0-D*S-M*X/Yxs;
    dy(2)=D*X0-D*X+M*X;
    dy(3)=D*P0-D*P+M*X*Ypx;
endfunction

S0=100;X0=0;P0=0;// g/l g/l  g/l
Ks=2;MI=0.4;Vol=300;F0=60;Yxs=0.5;Ypx=0.6;// g/l g/h h^-1  m3  m3/h

D=F0/Vol;

Xr=10;Sr=50;Pr=10;
y0=[Sr;Xr;Pr];x0=0;
t=100;
sol=ode(y0,x0,t,f);
disp(sol,"solucao")
