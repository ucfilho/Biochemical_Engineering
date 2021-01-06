clc
clear

function Derivada=f(t,Y,Num)

    
    for i=2:Num
      S=Y(i,1);X=Y(i,2);P=Y(i,3);
      Sb=Y(i-1,1);Xb=Y(i-1,2);Pb=Y(i-1,3);
      Derivada(i,1)=-U*(S-Sb)/deltaZ+X*MI*S/(S+Ks)*(1/Yxs)*(-1);
 	  Derivada(i,2)=-U*(X-Xb)/deltaZ+X*MI*S/(S+Ks);
 	  Derivada(i,3)=-U*(P-Pb)/deltaZ+X*MI*S/(S+Ks)*(Yps/Yxs)
    end
endfunction     


MI=0.4;Ks=2;Yps=0.4;Yxs=0.5;
Num=20;U=1;//m/s
S0=100;X0=10;P0=0;t0=0;Time=0.1:0.1:7;L=20000;

U=U*3600; // m/s
for i=1:Num
      Conc0(i,1)=S0;Conc0(i,2)=X0;Conc0(i,3)=P0;
end
deltaZ=L/Num
Ref=size(Time,2)
for k=1:Ref
    t=Time(k);
    Conc=ode(Conc0,t0,t,f);
    Y(k,:)=Conc(Num,:);
end
X=Y(:,2);S=Y(:,1);P=Y(:,3);

scf();a=get("current_axes");a.data_bounds=[0,0;10,100];plot(Time,S);xgrid(5);xtitle( 'Substrato' ) ;
xlabel("tempo(h)");ylabel("S(g/l)");
scf();a=get("current_axes");a.data_bounds=[0,0;10,100];plot(Time,X);xgrid(5);xtitle( 'Celulas' ) ;
xlabel("tempo(h)");ylabel("X(g/l)");
scf();a=get("current_axes");a.data_bounds=[0,0;10,100];plot(Time,P);xgrid(5);xtitle('Produto') ;
xlabel("tempo(h)");ylabel("P(g/l)");
