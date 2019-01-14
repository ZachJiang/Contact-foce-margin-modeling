%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Ritz method in one dimension with fmincon
%%%This file is created by Zachary Chunli JIANG
%%%Date: 2018/Dec/22
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


global steps;
global s L;

%s is in terms of the coordinate system of the object
%L is the total length of the paper strip
steps=80;
L=1; 
s=linspace(0,L,steps);


%%
options = optimoptions('fmincon', 'OutputFcn', @outfun,'Display','iter-detailed','Algorithm','interior-point');
options.MaxFunctionEvaluations = 1000000;
options.MaxIterations = 600;
options.OptimalityTolerance=1e-20;
options.ConstraintTolerance=1e-20;
options.StepTolerance=1e-40;

% Invoke optimization with a strating guess
a=[0.1,0.1,1,1,0.1,0.1,0.1,0.1,0.1,0.1];
[a,fval] = fmincon(@total_energy, a,[],[],[],[],[],[], @constraint_function, options);

figure
% plot(s,theta)
% hold on
%Set the basic function in terms of Fourier series
e1=1;
e2=s;
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

xc = cumtrapz(s, cos(theta));
yc = cumtrapz(s, sin(theta));

plot(xc, yc)


%%
function obj_value = total_energy(a)
%Set the basic function in terms of Fourier series
    global s L steps;            
    %Apply Fourier's form to approximate THETA function 
    e1=1;
    e2=s;
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

    dtheta_ds=diff(theta)./diff(s);
    
    %Currently only bending energy considered
    obj_value = trapz(s(1:end-1),(dtheta_ds).^2);

end
%%
function [c,c_eq] = constraint_function(a)
%Set the basic function in terms of Fourier series
    global s L;
    e1=1;
    e2=s;
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

    %end point tangent constraints
    c_eq_1 = theta(1);
    %c_eq_2 = theta(end)+pi/6;
    %end point constraints (x,y);
    x_e = 0.8*L;
    y_e = 0.2*L;
    c_eq_3 = trapz(s, cos(theta))-x_e;
    c_eq_4 = trapz(s, sin(theta))-y_e;

    yc = cumtrapz(s, sin(theta));
    
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
        
    visualize_plots(x)
    
    
end

%%
function visualize_plots(a)
%Set the basic function in terms of Fourier series

    global s L;
    e1=1;
    e2=s;
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

    xc = cumtrapz(s, cos(theta));
    yc = cumtrapz(s, sin(theta));
    
    plot(xc,yc)
    
    drawnow
    
    hold on
    
end
