clc
clear

function Derivada=f(t,Y)
      S=Y(1);X=Y(2);P=Y(3);
      Derivada(1,1)=1/U*X*MI*S/(S+Ks)*(1/Yxs)*(-1);
 	  Derivada(1,2)=1/U*X*MI*S/(S+Ks);
 	  Derivada(1,3)=1/U*X*MI*S/(S+Ks)*(Yps/Yxs)
endfunction     

U=1
MI=0.4;Ks=2;Yps=0.4;Yxs=0.5;
S0=100;X0=10;P0=0;z0=0;L=20000;
U=U*3600; // m/h
z=0.1:0.1:L;
Conc0(1)=S0;Conc0(2)=X0;Conc0(3)=P0;
Conc=ode(Conc0,z0,z,f);

Conc=Conc';
X=Conc(:,2);S=Conc(:,1);P=Conc(:,3);

scf();a=get("current_axes");a.data_bounds=[0,0;10,100];plot(z,S);xgrid(5);xtitle( 'Substrato' ) ;
xlabel("z(m)");ylabel("S(g/l)");
scf();a=get("current_axes");a.data_bounds=[0,0;10,100];plot(z,X);xgrid(5);xtitle( 'Celulas' ) ;
xlabel("z(m)");ylabel("X(g/l)");
scf();a=get("current_axes");a.data_bounds=[0,0;10,100];plot(z,P);xgrid(5);xtitle('Produto') ;
xlabel("z(m)");ylabel("P(g/l)");
