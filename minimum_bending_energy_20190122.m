%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Paper trajectory planning
%%%Ritz method in one dimension with fmincon
%%%This file is created by Zachary Chunli JIANG
%%%Date: 2018/Jan/22

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%Ref. One dimensional finite element Schaums ,F. Scheid,Numerical Analysis  
%Page 434
%In Hirai's IJRR Paper, basic function is expand in Fourier format
%fmincon is adopted to optimize the minimal energy

%%
%Clear the workspace parameters
clear all;


figure;

%%
global options;
options = optimoptions('fmincon', 'OutputFcn', @outfun,'Display','iter-detailed','Algorithm','interior-point');
options.MaxFunctionEvaluations = 1000000;
options.MaxIterations = 500;
options.OptimalityTolerance=1e-2;
options.ConstraintTolerance=1e-2;
options.StepTolerance=1e-8;

%%
% Invoke optimization with a strating guess
global steps;
global s L;
global x_end y_end;
global c3
global a a0 a1;
a=zeros(10,10);
% global concave_up_flex concave_down_flex total_flex_energy;
% global collect_x_end collect_y_end;
global kk;
%s is in terms of the coordinate system of the object
%L is the total length of the paper strip
steps=100;
L=120; 
s=linspace(0,L,steps);

%%
xo=0;
yo=0;
x0=(0:(L*1.2+160))-L*1.2-80;
y0=zeros(1,(L*1.2+161));
plot(x0,y0,'k');
axis([-L*1.2-80 80 -5 60]);
drawnow
hold on;
drawnow
scatter(xo,yo,'k','filled');
drawnow
hold on;
axis equal;
hold on;


kk=0:0.1:1;
c1=[0.12,1,1,1,1,0.5,0.4,1,0.21,1];
c2=[1,1,1,0.2,1,0.3,1,1,1,1];
c3=[0.2,0.2,0.9,0.9,0.1,0.1,0.1,0.1,0.1,0.1];
%%
%This is the main logic
scale=0.05;%% minimal 0.02

obj_value=zeros(1,10);
%start initial deformation optimization
x_end=L;
y_end=0;
a0=start_shape_optimization(c3);
a(1,:)=a0;
%start goal deformation optimization
x_end=0.9*L;
y_end=0.0*L;
a1=start_shape_optimization(c3);
a(10,:)=a1;
%start path planning

c=[c1;c2;c3];

c=start_trajectory_planning(c);


figure;
k=0;
for i=0:0.1:1
    k=k+1;
    a_i=(1-i)*a0+i*a1+c(1,:)*i*(1-i)+c(2,:)*i*i*(1-i)+c(3,:)*i^3*(1-i);
    theta_i=ritz_fourier_model(a_i);
    xc = L*cumtrapz(s/L, cos(theta_i))-L;
    yc = L*cumtrapz(s/L, sin(theta_i));
    plot(xc, yc);
    drawnow
    text(xc(end),yc(end),num2str(k));
    if k==4
        arc_length=trapz(xc,sqrt(1+(gradient(yc)./gradient(xc)).^2));
        disp(arc_length)
    end
    hold on;
end

grid on;
legend('1','2','3','4','5','6','7','8','9','10','11');
hold off;

%%
function new_a=start_trajectory_planning(c)
    global options;
 
    [c,fval] = fmincon(@min_max_energy, c,[],[],[],[],[],[], @planning_constrain_function, options);

    
    new_a=c;         

end


function obj=min_max_energy(c)
    global s L;
    global a0 a1 c3;
    global x_end y_end;
    obj_value=zeros(1,10);
    k=0;
    for i=0:0.1:1
        k=k+1;
        a_i=(1-i)*a0+i*a1+c(1,:)*i*(1-i)+c(2,:)*i*i*(1-i)+c(3,:)*i^3*(1-i);
        theta_i=ritz_fourier_model(a_i);
        
        
        xi = L*cumtrapz(s/L, cos(theta_i));
        yi = L*cumtrapz(s/L, sin(theta_i));
        x_end=xi(end);
        y_end=yi(end);

        a_i=start_shape_optimization(a_i);
        theta_i=ritz_fourier_model(a_i);
        dtheta_ds_i=diff(theta_i)./diff(s/L);
        theta_i=ritz_fourier_model(a_i);
        obj_value(k)= trapz(s(1:end-1)/L,(dtheta_ds_i).^2);    
        xc = L*cumtrapz(s/L, cos(theta_i))-L;
        yc = L*cumtrapz(s/L, sin(theta_i));
    
        plot(xc, yc);
        drawnow
        hold on;
    end
    obj=max(obj_value);
end

function [c,c_eq] = planning_constrain_function(a)
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
    xc = L*cumtrapz(s/L, cos(theta));
    
    c_ineq_1= -yc;
    c_ineq_2= -xc;
    c_eq = [c_eq_1];  
    c=[c_ineq_1;c_ineq_2];
end

function result_a=start_shape_optimization(az)

    global s L options;
    [az,fval] = fmincon(@total_energy, az,[],[],[],[],[],[], @constrain_function, options);
    theta=ritz_fourier_model(az);

    xc = L*cumtrapz(s/L, cos(theta))-L;
    yc = L*cumtrapz(s/L, sin(theta));

    plot(xc, yc);
    drawnow
    hold on;
    
    result_a=az;

end
function result_a=start_inter_shape_optimization(az)

    global s L options;
    [az,fval] = fmincon(@total_energy, az,[],[],[],[],[],[], @planning_constrain_function, options);
    theta=ritz_fourier_model(az);

    xc = L*cumtrapz(s/L, cos(theta))-L;
    yc = L*cumtrapz(s/L, sin(theta));

    plot(xc, yc);
    drawnow
    hold on;
    
    result_a=az;

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
    xc = L*cumtrapz(s/L, cos(theta));
    
    c_ineq_1= -yc;
    c_ineq_2= xc-x_e;
    c_eq = [c_eq_1; c_eq_3; c_eq_4];  
    c=[c_ineq_1;c_ineq_2];
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
     if 1-isempty(inflection_points)   
         distances = sqrt(sum(bsxfun(@minus, s/L, inflection_points(1)).^2,1));
         inflection_point1=find(distances==min(distances));
         flex_up=1/2*trapz(s(1:inflection_point1)/L,(dtheta_ds(1:inflection_point1)).^2);
     else 
         flex_up=-1;
     end
end

%%
function flex_down=concave_down_flex_energy(a)
    global s L;
    theta=ritz_fourier_model(a);
    dtheta_ds=diff(theta)./diff(s/L);
    inflection_points=find_cross(s(1:end-1)'/L,dtheta_ds',0);
    
    if 1-isempty(inflection_points)
        distances = sqrt(sum(bsxfun(@minus, s/L, inflection_points(1)).^2,1));
        inflection_point1=find(distances==min(distances));
        flex_down=1/2*trapz(s(inflection_point1:end-1)/L,(dtheta_ds(inflection_point1:end)).^2);
    else
        flex_down=-1;
    end
end

function k=curvature(a)
    global s L;
    theta=ritz_fourier_model(a);

    k=(diff(theta)./diff(s/L))/L;
end
