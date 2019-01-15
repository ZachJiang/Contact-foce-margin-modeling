%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Ritz method in one dimension with fmincon
%%%This file is created by Zachary Chunli JIANG
%%%Date: 2018/Jan/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%Ref. One dimensional finite element Schaums ,F. Scheid,Numerical Analysis  
%Page 434
%In Hirai's IJRR Paper, basic function is expand in Fourier format
%fmincon is adopted to optimize the minimal energy

%%
%Clear the workspace parameters
clear all;
clc;
close all;

figure;

%%
global options;
options = optimoptions('fmincon', 'OutputFcn', @outfun,'Display','iter-detailed','Algorithm','interior-point');
options.MaxFunctionEvaluations = 1000000;
options.MaxIterations = 600;
options.OptimalityTolerance=1e-20;
options.ConstraintTolerance=1e-20;
options.StepTolerance=1e-40;

%%
% Invoke optimization with a strating guess
global steps;
global s L;
global x_end y_end;
%s is in terms of the coordinate system of the object
%L is the total length of the paper strip
steps=100;
L=50; 
s=linspace(0,L,steps);

%%
xo=0;
yo=0;
x0=(0:(L*1.5+80))-L*1.5;
y0=zeros(1,(L*1.5+81));
plot(x0,y0,'k');
axis([-L*1.5 80 -20 60]);
drawnow
hold on;
drawnow
scatter(xo,yo,'k','filled');
drawnow
hold on;



%%
for theta1=0:1:0
    x1=30/12*theta1+50;
    z1=(theta1-100.6)/-0.8853;
    %disp('the current data set is: ');
    %dis=[x1,z1,theta1];
    %disp(dis);
    
    x=-x1;%0:10:-100
    z=z1+20; %90:1:110 and 20 is the offset
    theta=theta1*pi/180;%0:pi/180:12*pi/180/

    p=0.00:0.002:0.3; % pressure is in MPa
    gamma=45*pi/180;
    m=45;
    n=40;
    j=0:0.1:20;
    k=sqrt(20*20-(j-20).^2);
    l=72;
    sj=size(j,2);
    
    %point O1:
    x1=x;
    z1=z;
    delta=l.*(0.115*p+0.001177);
    phi=pi-delta-gamma+theta;
    r=1./(0.125*p+0.001177);
    %point O2:
    x2=x1+m*sin(theta)+n*cos(theta);
    z2=z1-m*cos(theta)+n*sin(theta);
    s1=cos(pi/2-gamma+theta);
    f1=sin(pi/2-gamma+theta);
    %point O3;
    x3=x2-r.*s1;
    z3=z2-r.*f1;
    %point E
    x4=x3+r.*cos(delta+gamma-pi/2-theta);
    z4=z3-r.*sin(delta+gamma-pi/2-theta);
    %point T
    x5=x4+j(sj).*sin(phi)-k(sj).*cos(phi);
    z5=z4-j(sj).*cos(phi)-k(sj).*sin(phi);
    for i1=1:size(p,2)
        scatter(x5(i1),z5(i1),3,'r');
        drawnow
        endpoint=[x5(i1),z5(i1)];
        
        if endpoint(1)<0 
            x_end=endpoint(1)+L;           
            if endpoint(2)<0
                y_end=0;
            else
                y_end= endpoint(2);
            end
    
            start_optimization();
        end
        
    end
    %legend('-DynamicLegend');
    
   
    hold all;
end
grid on;
hold off;
%%
function start_optimization()
    
    global s L options;
    a=[0.1,0.1,1,1,0.1,0.1,0.1,0.1,0.1,0.1];
    [a,fval] = fmincon(@total_energy, a,[],[],[],[],[],[], @constrain_function, options);
    theta=ritz_fourier_model(a);

    xc = L*cumtrapz(s/L, cos(theta))-L;
    yc = L*cumtrapz(s/L, sin(theta));

    plot(xc, yc);
    drawnow
    hold on;
    
    %calculate flex energies;
    concave_up_flex=concave_up_flex_energy(a);
    concave_down_flex=concave_down_flex_energy(a);
    total_flex_energy=concave_up_flex+concave_down_flex;
    disp('concave_up_energy:');
    disp(concave_up_flex);
    disp('concave_down_energy:');
    disp(concave_down_flex);
    disp('total_flex_energy:');
    disp(total_flex_energy);
end

%%
function obj_value = total_energy(a)
%Set the basic function in terms of Fourier series          
    %Apply Fourier's form to approximate THETA function 
    global s L;
    theta=ritz_fourier_model(a);

    dtheta_ds=diff(theta)./diff(s/L);
    
    %Currently only bending energy considered
    obj_value = trapz(s(1:end-1)/L,(dtheta_ds).^2);

end
%%
function [c,c_eq] = constrain_function(a)
%Set the basic function in terms of Fourier series
    global s L x_end y_end;
    theta=ritz_fourier_model(a);

    %end point tangent constraints
    c_eq_1 = theta(1);
    %c_eq_2 = theta(end)+pi/6;
    %end point constraints (x,y);
    x_e = x_end;
    y_e = y_end;
    c_eq_3 = L*trapz(s/L, cos(theta))-x_e;
    c_eq_4 = L*trapz(s/L, sin(theta))-y_e;

    yc = L*cumtrapz(s/L, sin(theta));
    
    c_ineq_1= yc;
    c_eq = [c_eq_1; c_eq_3; c_eq_4];  
    c=[-c_ineq_1];
end

%%
function stop = outfun(x, optimValues, state)
    stop = false;

    global history;
    if isequal(state,'iter')
          history = [history; x];
    end
        
    %visualize_plots(x)
    
    
end
%%
function visualize_plots(a)
%Set the basic function in terms of Fourier series
    
    global s L;
    theta=ritz_fourier_model(a);

    xc = L*cumtrapz(s/L, cos(theta));
    yc = L*cumtrapz(s/L, sin(theta));
    
    plot(xc,yc)
    
    drawnow
    
    hold on
    
end

%%
function theta=ritz_fourier_model(a)

    global s L;
    e1=1;
    e2=s/L;
    e3=sin(2*pi.*s/L);
    e4=cos(2*pi.*s/L);
    e5=sin(4*pi.*s/L);
    e6=cos(4*pi.*s/L);
    e7=sin(6*pi.*s/L);
    e8=cos(6*pi.*s/L);
    e9=sin(8*pi.*s/L);
    e10=cos(8*pi.*s/L);
    
    theta=a(1)*e1+a(2)*e2+a(3)*e3+a(4)*e4+a(5)*e5 ...
    +a(6)*e6+a(7)*e7+a(8)*e8+a(9)*e9+a(10)*e10;

end

%%
function theta=symbolic_theta(a)
    
    global L;
    syms t;
    e1=1;
    e2=t/L;
    e3=sin(2*pi.*t/L);
    e4=cos(2*pi.*t/L);
    e5=sin(4*pi.*t/L);
    e6=cos(4*pi.*t/L);
    e7=sin(6*pi.*t/L);
    e8=cos(6*pi.*t/L);
    e9=sin(8*pi.*t/L);
    e10=cos(8*pi.*t/L);
    
    r(t)=a(1)*e1+a(2)*e2+a(3)*e3+a(4)*e4+a(5)*e5 ...
    +a(6)*e6+a(7)*e7+a(8)*e8+a(9)*e9+a(10)*e10;
    
    theta=matlabFunction(r);
    
end


%%
function flex_up=concave_up_flex_energy(a)
%tutorial:https://www.mathworks.com/help/symbolic/examples/maxima-minima-and-inflection-points.html
    global s L;
    theta=ritz_fourier_model(a);
    dtheta_ds=diff(theta)./diff(s/L);
    inflection_points=find_cross(s(1:end-1)'/L,dtheta_ds',0);
    
    distances = sqrt(sum(bsxfun(@minus, s/L, inflection_points(1)).^2,1));
    inflection_point1=find(distances==min(distances));
    
    flex_up=1/2*trapz(s(1:inflection_point1)/L,(dtheta_ds(1:inflection_point1)).^2);
end

function flex_down=concave_down_flex_energy(a)
    global s L;
    theta=ritz_fourier_model(a);
    dtheta_ds=diff(theta)./diff(s/L);
    inflection_points=find_cross(s(1:end-1)'/L,dtheta_ds',0);
    
    distances = sqrt(sum(bsxfun(@minus, s/L, inflection_points(1)).^2,1));
    inflection_point1=find(distances==min(distances));
    
    flex_down=1/2*trapz(s(inflection_point1:end-1)/L,(dtheta_ds(inflection_point1:end)).^2);
end
