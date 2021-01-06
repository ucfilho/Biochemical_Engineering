clc;clear;
U=0.6; // m/h
flow=30; // m3/h
theta=5;// h
Area=flow/U;disp(Area,"Area (m2)=")
D=(4*Area/%pi)^0.5;disp(D,"Diametro(m)=")
Vol=flow*theta;disp(Vol,"Volume (m3)=")
Altura=Vol/Area;disp(Altura,"Altura(m)=")
