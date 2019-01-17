%Here we compute the minimum bending energy curve
%with fixed end-point constraints. We use nonlinear
%optimization to minimize the bending energy.

%Author: Abdullah Nazir

clear all

global var_s

intervals = 70; %number of sample point intervals

arc_length = 1; %Fixed curve length constraint

var_s = linspace(0, arc_length, intervals); %arclength variable

var_theta_init = 1*ones(1, length(var_s)) +...
                1*sin(2*pi*var_s/var_s(end)) +...
                1*cos(2*pi*var_s/var_s(end)) +...
                1*sin(4*pi*var_s/var_s(end)) +...
                1*cos(4*pi*var_s/var_s(end)); %Initial guess for theta(s)



options = optimoptions('fmincon',...
                    'OutputFcn', @outfun,...
                    'Display',...
                    'iter-detailed',...
                    'Algorithm',...
                    'interior-point');
                
options.MaxFunctionEvaluations = 100000000;
% options.MaxIterations = 200;
options.OptimalityTolerance= 1e-2; % important parameter
% options.ConstraintTolerance= 1e-10;
options.StepTolerance = 1e-3; %important parameter. try commenting this and
%change optimalityTolerance to 1e-2



var_theta = fmincon(@objectiveFunction,...
                var_theta_init,...
                [],[],[],[],[],[],...
                @constraintFunctions,...
                options);



figure
[xc, yc] = arcLengthToCartesian(var_theta, var_s);

plot(xc, yc, 'Color', 'r')
hold on
% plot(var_s, var_theta, 'Color', 'g')
% hold on
% plot(var_s, gradient(var_theta)./gradient(var_s), 'Color', 'b')



U_flex = computeFlexuralEnergy(var_theta, var_s)

%Compute inflextion point and inflexion index of the array
[inflexion_point, inf_idx] = computeInflexionPoint(xc, yc);

U_flex_CU = computeFlexuralEnergy(var_theta(1:inf_idx), var_s(1:inf_idx))

U_flex_CD = computeFlexuralEnergy(var_theta(inf_idx+1:end), var_s(inf_idx+1:end))

% 
% dydx = gradient(yc)./gradient(xc);
% d2ydx = gradient(dydx)./gradient(xc);

function fvalue = objectiveFunction(var_theta)

    global var_s

    [dthetads, d2thetads] = computeDifferentials(var_theta, var_s);
    
    fvalue = trapz(var_s, dthetads.^2);
    
end


function [c, c_eq] = constraintFunctions(var_theta)
    %This function defines equality and inequality constraints on the 
    %optimization problem. 
    
    global var_s
    
    % end point tangent constraint
    c_eq_1 = var_theta(1);
%     c_eq_2 = var_theta(end);

    x_end = 0.9;
    y_end = 0.05;
    
    % end point position constraint (arc length coordinates)
    c_eq_3 = trapz(var_s, cos(var_theta))-x_end;
    c_eq_4 = trapz(var_s, sin(var_theta))-y_end;
    
    c_eq = [c_eq_1; c_eq_3; c_eq_4];
    
    c = [];

end



function [dthetads, d2thetads] = computeDifferentials(var_theta, var_s)

    dthetads = gradient(var_theta)./gradient(var_s);
    
    d2thetads = gradient(dthetads)./gradient(var_s);

end


function [xc, yc] = arcLengthToCartesian(var_theta, var_s)

    %This functions converts arc length coordinates to cartesian
    %coordinates of the curve

    xc = cumtrapz(var_s, cos(var_theta));
    yc = cumtrapz(var_s, sin(var_theta));
end

function U_flex = computeFlexuralEnergy(var_theta, var_s)
    %this function computes flexural energy of a curve based on the 
    %curvature, without factoring the rigidity constraint.
    
    [dthetads, d2thetads] = computeDifferentials(var_theta, var_s);
    
    U_flex = trapz(var_s, dthetads.^2);
    

end

function [inflexion_point, inflexion_idx] = computeInflexionPoint(xc, yc)

    dydx = gradient(yc)./gradient(xc);
    d2ydx = gradient(dydx)./gradient(xc);   
    
    %Too simplistic calculation for inflexion point. Account for multiple
    %inflexion points.
    id = sign(d2ydx);
    inflexion_idx = strfind(id,[1 -1]) + 1;
    inflexion_point = xc(inflexion_idx);

    
end





function stop = outfun(var_theta, optimValues, state)
    stop = false;

    global history;
    if isequal(state,'iter')
          history = [history; var_theta];
    end
        
    visualize_plots(var_theta)
    
end

function visualize_plots(var_theta)

    global var_s
    
    plot(var_s, var_theta)
    drawnow
    hold on
end
