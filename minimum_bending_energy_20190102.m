% currently there some bugs with this code 2019_1_4
% 1. the length of the paper strip is not fixed
% 2. the object function is not based on real parameters
% 3. should add inequality constrains for some cases to cover

clear all
clc 
global t history steps
global arc_length

arc_length = 3;
steps = 100;
history = zeros(1,steps);

t_min = 0;
t_max = 1;

t = linspace(t_min,t_max,steps);

% 
 x0 = 1*ones(1,steps);
 x0(1)=0;
 x0(end)=0;

 x0 = (0.5-abs(t_max/2-t))*arc_length;

 %x0 = 1-cos(t*2*pi/t_max);

options = optimoptions('fmincon', 'OutputFcn', @outfun,'Display','iter-detailed','Algorithm','sqp');
options.MaxFunctionEvaluations = 10000000;
options.MaxIterations = 2000;
options.StepTolerance = 1e-40;
options.OptimalityTolerance=1e-40;
options.ConstraintTolerance=1e-40;

x = fmincon(@bending_energy, x0,[],[],[],[],[],[], @constraint_function, options);


[dxdt, d2xdt] = compute_derivatives(x,t);

figure
plot(t,x)
% title('t,x relation')
% hold on
% plot(t, dxdt)
% title('t, dxdt relation')
% hold on
% plot(t, d2xdt)
% title('t,d2xtx relation')

% axis([0 3 0 2])


function obj_value = bending_energy(x)

    global t history 
    
    
    [dxdt, d2xdt] = compute_derivatives(x,t);    
    
    del_x = x - history(end,:);
       
        
    obj_value = trapz(t,d2xdt.^2);


end

function [c,c_eq] = constraint_function(x)

    
    global t arc_length
    
    [dxdt, d2xdt] = compute_derivatives(x,t);    

    
    c_eq_1 = x(1);
    c_eq_2 = x(end);

    c_eq_3 = dxdt(1);
    c_eq_4 = dxdt(end)+pi/2;
    
    int_ds = trapz(t, sqrt(1+dxdt.^2));
    c_eq_5 = int_ds-arc_length;
    
    
    c_eq = [c_eq_1; c_eq_2; c_eq_3; c_eq_4; c_eq_5];

    c=[];
end

function [dxdt, d2xdt] = compute_derivatives(x,t)

    dt = gradient(t);
    dx = gradient(x);
    
    dxdt = dx./dt;
    
    d2xdt = gradient(dxdt)./dt;
    
end





function stop = outfun(x, optimValues, state)
    stop = false;

    global history;
    if isequal(state,'iter')
          history = [history; x];
    end
        
    visualize_plots(x)
    
    
end

function visualize_plots(x)

    global t
    
    plot(t, x)
    
    drawnow
    
    hold on
    
end

