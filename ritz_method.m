clear all

%optimization variable 'theta_coeff' is a collection of function coefficients.
%each function linearly independent

global s

% set the step size of each iteraction
steps = 500;

s_min = 0; 
s_max = 1;

s = linspace(s_min,s_max,steps);

theta_coeff0 = ones(1,5);

% the optimization settings
options = optimoptions('fmincon', 'OutputFcn', @outfun,'Display','iter-detailed','Algorithm','sqp');
options.MaxFunctionEvaluations = 10000000;
options.MaxIterations =  30;
options.StepTolerance=1e-20;
options.OptimalityTolerance=1e-10;
options.ConstraintTolerance=1e-10;

theta_coeff = fmincon(@all_kinds_energy, theta_coeff0,[],[],[],[],[],[], @constraint_function, options);

theta = coeffToFunction(theta_coeff, s);

[xc, yc] = toCartesian(theta, s);

figure
daspect([1 1 1])
plot(xc, yc)

% objective function setting
function obj_value = all_kinds_energy(theta_coeff)

    global s
    
    theta = coeffToFunction(theta_coeff, s);
    
    dthetads = compute_derivatives(theta,s);
    % bending energy          
    bending_value = trapz(s,(dthetads).^2);
    
    [xt, yt] = toCartesian(theta, s);
    % gravity energy
    grav_value=trapz(s,yt);
    
    obj_value=grav_value+bending_value;
end

% constrain functions 
% c_eq is the equality constrain set
% c is the inequality constrain set that c<=0
function [c,c_eq] = constraint_function(theta_coeff)

    global s
    
    theta = coeffToFunction(theta_coeff, s);
    
    %%%%%%%%%%%%%%%%%%%%%end point tangent constraints%%%%%%%%%%%%%%%%%%%%
    theta_s = 0;  
    theta_e = -60*pi/180;   
    
    c_eq_1 = theta(1)-theta_s;
    c_eq_2 = theta(end)-theta_e;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%%%%%%%%%%%%%%%%%%%end point constraints (x,y)%%%%%%%%%%%%%%%%%%%%%%
    x_e = 0.8;
    y_e = 0;

    c_eq_3 = trapz(s, cos(theta))-x_e;
    c_eq_4 = trapz(s, sin(theta))-y_e;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    [xc, yc] = toCartesian(theta, s);
    
    c_eq_5 = trapz(xc, sqrt(1+(gradient(yc)./gradient(xc)).^2))-s(end);
        
    c_eq = [c_eq_1; c_eq_2; c_eq_3;c_eq_4];

    c=[-yc];
end

% expand thetha in 
function theta = coeffToFunction(theta_coeff, s)
    
    theta = theta_coeff(1)*ones(1, length(s)) ...
            + theta_coeff(2)*sin(2*pi*s/s(end)) ...
            + theta_coeff(3)*cos(2*pi*s/s(end));
%             + theta_coeff(4)*sin(4*pi*s/s(end)) ...
%             + theta_coeff(5)*cos(4*pi*s/s(end));
            
end

% derivative of theta
function dthetads = compute_derivatives(theta,s)


    dthetads = gradient(theta)./gradient(s);
    
    
end

% transform from local (theta,s) frame alongs the shape of curve to global cartesian frame
function [xc, yc] = toCartesian(theta, s)

    xc = cumtrapz(s, cos(theta));
    yc = cumtrapz(s, sin(theta));
end

% stopping criteria for iterations
function stop = outfun(theta_coeff, optimValues, state)
    stop = false;

    global history;
    if isequal(state,'iter')
          history = [history; theta_coeff];
    end
        
    visualize_plots(theta_coeff)
end

%plot out the results
function visualize_plots(theta_coeff)

    global s
    
    theta = coeffToFunction(theta_coeff, s);
    
    plot(s, theta)
    
    drawnow
    
    hold on
    
end

