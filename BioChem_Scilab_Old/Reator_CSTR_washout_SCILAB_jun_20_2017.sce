clc
clear

Ks=2;MI=0.4;Vol=300;F0=60;Yxs=0.5;Ypx=0.6;// g/l g/h h^-1  m3  m3/h
S0=100;X0=0;P0=0;// g/l g/l  g/l

function dy=f(x,y)
    S=y(1);X=y(2);P=y(3);
    M=MI*S/(S+Ks)
    dy(1)=D*S0-D*S-M*X/Yxs;
    dy(2)=D*X0-D*X+M*X;
    dy(3)=D*P0-D*P+M*X*Ypx;
endfunction

Num=100;
for i=1:Num
    D=0.1+0.6*(i-1)/(Num-1)
    Di(i)=D;
    Xr=10;Sr=50;Pr=10;
    y0=[Sr;Xr;Pr];x0=0;
    t=500;
    sol=ode(y0,x0,t,f);
    Si(i)=sol(1)
end





scf()
plot(Di,Si);xgrid(5);xtitle( 'Substrato CSTR' ) ;

xlabel("$ D (h^{-1})$");ylabel("S(g/l)");

a=get("current_axes");
a.font_size=3; 
x_label=a.x_label;
x_label.font_size= 5;
a.children.children.thickness=4;
y_label=a.y_label;
y_label.font_size= 5;
t=a.title;
t.foreground=9;
t.font_size=4;
t.font_style=5;

