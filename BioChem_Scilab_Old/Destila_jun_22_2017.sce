clc;clear;clf

//============================================================================

//  PARTE DE ENTRADA DE DADOS

//============================================================================

//equilibrium data-ethanol
xeq_d =[0.0010 0.0061 0.0145 0.0237 0.0310 0.0490 0.0652 0.0968 0.1394 0.3261 0.4635 0.5413 0.6856 0.7760 0.8403 0.9037 0.9725 0.9804];

yeq_d =[0.0047 0.0721 0.1539 0.2301 0.2851 0.3559 0.4181 0.4534 0.5314 0.6047 0.6518 0.6751 0.7451 0.8005 0.8457 0.9010 0.9721 0.9774];


xd = 0.84; //  xD     - distillate mole fraction
xf = 0.0400722;//  xF     - feed mole fraction
xw=0.0039370;//  xW     - bottoms mole fraction
q  = 1;//  q      - mole fraction of liquid in the feed
RR =3; //  R      - reflux ratio

AJUSTE=3; 
// OBS: AJUSTE=1 SPLINLINE, AJUSTE=2 y=((a+b*x)./(1+c*x+d*x.^2)), 
//      AJUSTE=3 y=((a1+b1*x+c1*x.^2)./(1+d1*x+e1*x.^2+f1*x.^3))
//      AJUSTE=4 raw data


//============================================================================
//  PARTE DE CALCULO DO CODIGO :
//  Programa adaptados de codigo de autoria do  Jakub Kopáč
//============================================================================

info_level = 1; // 0 = no info printed,  1 = calculated minimal reflux printed  2 = fsolve informations printed 
a = 0; b = 1;        // interval of interpolation
m = 500;             // discretisation for evaluation
xeq= linspace(a,b,m);

printf(" \n ---Inicio do programa--------- \n")
printf(" \n \n McCabe-Thiele method \n \n")
printf(" \n xd (distillate mole fraction) = %g \n", xd);
printf(" \n xf (feed mole fraction) = %g \n",xf);
printf(" \n xw ( bottoms mole fraction) = %g \n",xw);
printf(" \n q ( mole fraction of liquid in the feed) = %g \n",q);
printf(" \n RR ( reflux ratio) = %g \n",RR);


function y=data_fit(x,coe)
    if AJUSTE==2 then
        a=coe(1);b=coe(2);c=coe(3);d=coe(4);y=((a+b*x)./(1+c*x+d*x.^2));
    elseif AJUSTE==3 then
        a1=coe(1);b1=coe(2);c1=coe(3);d1=coe(4);e1=coe(5);f1=coe(6);
        y=((a1+b1*x+c1*x.^2)./(1+d1*x+e1*x.^2+f1*x.^3))
    end
endfunction

function Erro=Model_Error(coe,z)
    x=z(1);y=z(2);Erro=y-data_fit(x,coe);
endfunction



if AJUSTE==1 then
    yeq = interp(xeq, xeq_d,yeq_d, splin(xeq_d,yeq_d,"not_a_knot"));
elseif AJUSTE==2 then
    Ainit=[1;1;1;1];z=[xeq_d;yeq_d];[popt,err]=datafit(Model_Error,z,Ainit);
    for i=1:m yeq(i)=data_fit(xeq(i),popt); end
elseif AJUSTE==3 then
    Ainit=[1;1;1;1;1;1];z=[xeq_d;yeq_d];[popt,err]=datafit(Model_Error,z,Ainit);
    for i=1:m yeq(i)=data_fit(xeq(i),popt); end
elseif AJUSTE==4 then
    xeq=xeq_d;yeq=yeq_d;
end

// 

function mct_main(xD,xF,xW,q,R,VLE,isInfo)
    // ### DEFINITION OF SUBFUNCTIONS ###
    function abort_for_error(str)
        disp(sprintf("McCabe-Thiele error: %s",str))
        disp("# PROGRAM ABORTED #")
        abort
    endfunction
    function y = sub_equilibrium(x,xy_eq_data)
        //for any x in range <0,1> interpolate the coresponding y
        //if or([x<0,x>1]) then
        //    abort_for_error("x value out of range <0,1>")
        //end
        y = interpln(xy_eq_data,x)
    endfunction
    function xiyi = sub_minreflux (xy)
        //find minimal reflux as 2 NAEs for intersection point of q-line and equilibrium line
        x   = xy(1);
        y   = xy(2);
        yeq     = sub_equilibrium(x,VLE);
        xiyi(1) = (y*(q-1)-(q*x-xF)); //q-line
        xiyi(2) = y-yeq;              //equlibrium line
    endfunction
    function xiyi = sub_intersection(xy)
        x   = xy(1);
        y   = xy(2);
        xiyi(1) = y-(x*R/(R+1)+xD/(R+1)) //operating line for striping section
        xiyi(2) = y*(q-1)+xF-x*q         //q-line
    endfunction
    function fx = sub_xstep(x,y)
        //hladam x-ovu suradnicu pre vodorovnu ciaru
        fx = y-interpln(VLE,x)
    endfunction
    function fy = sub_ystep(y,x)
        //hladam y-ovu suradnicu pre zvilsu ciaru
        fy = y-interpln([x_work_lines;y_work_lines],x)
    endfunction
    function sub_graphics(xF,xD,xW,xeq,yeq,xwl,ywl,steps,n_steps)
        label_str_x = "x = mole fraction of lighter component in liquid phase";
        label_str_y = "y = mole fraction of lighter component in vapour phase";

        //q-line
        plot([xF xyRmin(1)],[xF xyRmin(2)])
        //working lines
        plot(xwl,ywl,"r")
        //diagonal + equilibrium line
        plot([xW xW],[0 xW],"k",[xF xF],[0 xF],"k",[xD xD],[0 xD],"k")
        plot([0 1],[0 1],"k",xeq,yeq,"b")
        //names of axis+title
        xlabel(label_str_x,"fontsize",2)
        ylabel(label_str_y,"fontsize",2)
        //title(label_title,"fontsize",3)

        //names of axis+title
        title(sprintf("McCabe-Thiele method, number of theoretical plates = %g",n_steps),"fontsize",3)
        printf("\n \n McCabe-Thiele method, number of theoretical plates = %g",n_steps);
        plot(steps(:,1),steps(:,2),"m")
    endfunction
    function x = sub_bisekcia(fun,funpar)
        a=0; b=1; tol=1e-8;
        fa = fun(a,funpar);
        fb = fun(b,funpar);

        if (fa*fb)>0 then
            abort_for_error("this should never happen...")
        end
        while (b-a) > tol
            x  = (a+b)/2;
            fx = fun(x,funpar);
            if fx==0 then
                return
            end
            if fa*fx<0 then
                b  = x;
            else
                a  = x;
                fa = fx;
            end
        end
    endfunction
    function [stepspoints,steps] = sub_stepskernel()
        steps_limit = 100
        x(1) = xD
        y(1) = xD
        ii    = 1
        steps = 0
        while x($)>xW
            ii    = ii+1
            x(ii) = sub_bisekcia(sub_xstep,y(ii-1))
            y(ii) = y(ii-1)

            ii    = ii+1
            x(ii) = x(ii-1)
            y(ii) = sub_bisekcia(sub_ystep,x(ii-1))
            steps = steps+1

            if steps>steps_limit then
                mprintf("\nWARNING: Current number of plates is %g",steps)
                WhatToDo = input("Do you want to continue ? (0/1)..")
                if WhatToDo==0 then
                    mprintf("\nCalculation CANCELLED - actual results will be ploted")
                    stepspoints = [x,y]
                    return
                elseif WhatToDo==1
                    steps_limit = steps_limit+30
                end
            end
        end

        y($) = 0    
        stepspoints = [x,y]
    endfunction

    if or([xd<xf,xd<xw,xf<xw]) then
        abort_for_error("Incorrect definition of mole fractions")
    end
    if interpln(VLE,xD)<xD then
        abort_for_error("Equilibrium y for xD is lower than xD")
    end

    //find minimal reflux ratio intersection points
    [xyRmin,fx1,info1] = fsolve([xF;xF],sub_minreflux);
    if isInfo>1 then
        mprintf("\nsub_minreflux")
        mprintf("\n   fx(1) = %1.2e\n   fx(2) = %1.2e",fx1(1),fx1(2))
        mprintf("\n   fsolve termination indicator = %d\n",info1)
    end

    //calculate minimal reflux from intersection point (I)
    // R_min = (xD-yI)/(yI-xI)
    R_min = (xD-xyRmin(2))/(xyRmin(2)-xyRmin(1));
    if isInfo>0 then
        mprintf("\n \n Calculated minimal reflux ratio = %1.2f",R_min)
    end
    if R<=R_min then
        abort_for_error("Reflux lower than calculated minimal reflux")
    end
    //working lines intersection
    [xyInter,fx2,info2] = fsolve([xF;xF],sub_intersection)
    if isInfo>1 then
        mprintf("\n\nsub_intersection")
        mprintf("\n   fx(1) = %1.2e\n   fx(2) = %1.2e",fx2(1),fx2(2))
        mprintf("\n   fsolve termination indicator = %d",info2)
    end
    //define x-y coordinates for working lines
    x_work_lines = [0 xW xyInter(1) xD 1];
    y_work_lines = [0 xW xyInter(2) xD 1];

    [stepspoints,n_steps] = sub_stepskernel();

    sub_graphics(xF,xD,xW,VLE(1,:),VLE(2,:),x_work_lines,y_work_lines,stepspoints,n_steps)
endfunction

xeq=xeq(:)';yeq=yeq(:)';
VLE = [xeq;yeq]; //  VLE    - vapour-liquid data
mct_main(xd,xf,xw,q,RR,VLE,info_level)

plot(xeq_d,yeq_d, 'o');
printf("\n \n ---Fim do programa---------\n")
