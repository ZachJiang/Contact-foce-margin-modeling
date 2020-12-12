%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Ritz method in one dimension with fmincon
%%%This file is created by Zachary Chunli JIANG
%%%Date: 2018/Jan/14
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

%%
options = optimoptions('fmincon', 'OutputFcn', @outfun,'Display','iter-detailed','Algorithm','interior-point');
options.MaxFunctionEvaluations = 1000000;
options.MaxIterations = 1000000;
options.OptimalityTolerance=1e-20;
options.ConstraintTolerance=1e-20;
options.StepTolerance=1e-40;

% Invoke optimization with a strating guess
global steps;
global s L;

%s is in terms of the coordinate system of the object
%L is the total length of the paper strip
steps=100;
L=1; 
s=linspace(0,L,steps);

a=[0.1,0.1,1,1,0.1,0.1,0.1,0.1,0.1,0.1];
[a,fval] = fmincon(@total_energy, a,[],[],[],[],[],[], @constraint_function, options);


theta=ritz_fourier_model(a);

xc = cumtrapz(s, cos(theta));
yc = cumtrapz(s, sin(theta));

figure
plot(xc, yc);
%calculate flex energies;

concave_up_flex=concave_up_flex_energy(a);
concave_down_flex=concave_down_flex_energy(a);
total_flex_energy=concave_up_flex+concave_down_flex;
disp('concave_up_energy:');
disp(concave_up_flex);
disp('concave_down_energy:');
disp(concave_down_flex);
disp('total_flex_energy:');
disp(total_flex_energy)



%%
function obj_value = total_energy(a)
%Set the basic function in terms of Fourier series          
    %Apply Fourier's form to approximate THETA function 
    global s;
    theta=ritz_fourier_model(a);

    dtheta_ds=diff(theta)./diff(s);
    
    %Currently only bending energy considered
    obj_value = trapz(s(1:end-1),(dtheta_ds).^2);

end
%%
function [c,c_eq] = constraint_function(a)
%Set the basic function in terms of Fourier series
    global s L;
    theta=ritz_fourier_model(a);

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
    
    global s;
    theta=ritz_fourier_model(a);

    xc = cumtrapz(s, cos(theta));
    yc = cumtrapz(s, sin(theta));
    
    plot(xc,yc)
    
    drawnow
    
    hold on
    
end

%%
function theta=ritz_fourier_model(a)

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

end

%%
function theta=symbolic_theta(a)
    
    global L;
    syms t;
    e1=1;
    e2=t;
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
    global s;
    theta=ritz_fourier_model(a);
    dtheta_ds=diff(theta)./diff(s);
    inflection_points=find_cross(s(1:end-1)',dtheta_ds',0);
    
    distances = sqrt(sum(bsxfun(@minus, s, inflection_points(1)).^2,1));
    inflection_point1=find(distances==min(distances));
    
    flex_up=1/2*trapz(s(1:inflection_point1),(dtheta_ds(1:inflection_point1)).^2);
end

function flex_down=concave_down_flex_energy(a)
    global s;
    theta=ritz_fourier_model(a);
    dtheta_ds=diff(theta)./diff(s);
    inflection_points=find_cross(s(1:end-1)',dtheta_ds',0);
    
    distances = sqrt(sum(bsxfun(@minus, s, inflection_points(1)).^2,1));
    inflection_point1=find(distances==min(distances));
    
    flex_down=1/2*trapz(s(inflection_point1:end-1),(dtheta_ds(inflection_point1:end)).^2);
end

