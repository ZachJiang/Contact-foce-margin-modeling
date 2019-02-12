%Here we compute the minimum bending energy curve
%with fixed end-point constraints. We use nonlinear
%optimization to minimize the bending energy.

clear all

global var_s x_end y_end

paper_length = 0.125;
paper_width = 0.050;
flex_rigidity_per_lw = 0.01615; %Unit Newtons
flex_rigidity = flex_rigidity_per_lw*paper_length*paper_width

return 

% flex_rigidity = 1;

intervals = 70; %number of sample point intervals

arc_length = paper_length; %Fixed curve length constraint

var_s = linspace(0, arc_length, intervals); %arclength variable

var_theta_init = 1*ones(1, length(var_s)) +...
                1*sin(2*pi*var_s/var_s(end)) +...
                1*cos(2*pi*var_s/var_s(end)) +...
                1*sin(4*pi*var_s/var_s(end)) +...
                1*cos(4*pi*var_s/var_s(end)); %Initial guess for theta(s)


options = optimoptions('fmincon',...
                    'OutputFcn', @outfun,...
                    'Algorithm',...
                    'interior-point');
                
options.MaxFunctionEvaluations = 100000000;
% options.MaxIterations = 200;
options.OptimalityTolerance= 1e-2; % important parameter
% options.ConstraintTolerance= 1e-10;
options.StepTolerance = 1e-3; %important parameter. try commenting this and
%change optimalityTolerance to 1e-2

x_bound = [0.3*arc_length, 0.98*arc_length];
y_bound = [0.0,0.45*arc_length];

% sizey = 20;
% sizex = ceil(sizey*((x_bound(2)-x_bound(1))/(y_bound(2)-y_bound(1))));


size = 20;

[X, Y] = sampleEnergyDomain(x_bound, y_bound, size);

U_flex_matrix = zeros(size);
U_flex_CU_matrix = zeros(size);
U_flex_CD_matrix = zeros(size);

CoF_matrix = zeros(size);


ones_lt_flipped = flip(tril(ones(size),1));
ut_idxes = find(ones_lt_flipped>0);

 
figure
axis equal
% axis([0 1 0 1])
% plot_boundary(arc_length)

for xidx = ut_idxes.'

        x_end = X(xidx);
        y_end = Y(xidx);
        

        [var_theta,fval,exitflag,output,lambda,grad,hessian] = fmincon(@objectiveFunction,...
                        var_theta_init,...
                        [],[],[],[],[],[],...
                        @constraintFunctions,...
                        options);

        [xc, yc] = arcLengthToCartesian(var_theta, var_s);
        plot(1000*(xc-paper_length),1000*yc)
        drawnow
        hold on

        U_flex_matrix(xidx) = flex_rigidity*computeFlexuralEnergy(var_theta, var_s);
        [U_flex_CU_matrix(xidx), U_flex_CD_matrix(xidx)] = compute_CUP_CAP_energies(var_theta, var_s);


        CoF_matrix(xidx) = computeCoF(lambda.eqnonlin(2),lambda.eqnonlin(3), var_theta(end));


%             %put force and normal vectors inside a function.
%             constraint_f_vec = [lambda.eqnonlin(2),lambda.eqnonlin(3)];
%             constraint_f_vec = 0.05*constraint_f_vec/norm(constraint_f_vec);
% 
%             quiver(x_end,y_end,constraint_f_vec(1),constraint_f_vec(2), 'Color', 'm')
%             drawnow
%             hold on
% 
%             if var_theta(end)<=0
%                 normal_vec = 0.05*[-sin(var_theta(end)), cos(var_theta(end)), 0];
%             else
%                 normal_vec = 0.05*[sin(var_theta(end)), -cos(var_theta(end)), 0];
%             end
% 
%             quiver(x_end,y_end,normal_vec(1),normal_vec(2), 'Color', 'k')
%             drawnow
%             hold on



end



%Convert zeros to NaN
U_flex_matrix(U_flex_matrix==0) = NaN; 
U_flex_CU_matrix(U_flex_CU_matrix==0) = NaN; 
U_flex_CD_matrix(U_flex_CD_matrix==0) = NaN; 

CoF_matrix(CoF_matrix==0) = NaN;


% surf(1000*(X-paper_length),1000*Y,U_flex_matrix)
% hold on
surf(1000*(X-paper_length),1000*Y,U_flex_matrix)
hold on


hx = (x_bound(2)-x_bound(1))/(size-1);
hy = (y_bound(2)-y_bound(1))/(size-1);

[gx,gy] = computeNumericalGradient(U_flex_matrix, hx, hy);

figure
axis equal

%start points on upper triangular line
% startx = linspace(0.40, 0.90, size);
% starty = -1*startx + 0.90; 
% streamline(X,Y,gx,gy,startx,starty)
% hold on
quiver(X,Y,-gx,-gy)

save('bending_curve_lambda_5', 'X', 'Y',...
    'U_flex_matrix', 'U_flex_CU_matrix', 'U_flex_CD_matrix',...
    'CoF_matrix','hx','hy')

% save('energy_function', 'X', 'Y', 'U_flex_matrix')

function fvalue = objectiveFunction(var_theta)

    global var_s

    [dthetads, d2thetads] = computeDifferentials(var_theta, var_s);
    
    fvalue = trapz(var_s, dthetads.^2);
    
end


function [c, c_eq] = constraintFunctions(var_theta)
    %This function defines equality and inequality constraints on the 
    %optimization problem. 
    
    global var_s x_end y_end
    
    % end point tangent constraint
    c_eq_1 = var_theta(1);
%     c_eq_2 = var_theta(end);

    % end point position constraint (cartesian)
%     x_end = 0.7;
%     y_end = 0.7;
    
    % end point position constraint (arc length coordinates)
    c_eq_3 = trapz(var_s, cos(var_theta))-x_end;
    c_eq_4 = trapz(var_s, sin(var_theta))-y_end;
    
    c_eq = [c_eq_1; c_eq_3; c_eq_4];
    
    c = [];

end

function plot_boundary(arc_length)

    t = linspace(0,2*pi/3,100);
    xt = arc_length*cos(t);
    yt = arc_length*sin(t);
    plot(xt, yt)
    hold on

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

function [gx,gy] = computeNumericalGradient(F, hx, hy)

    %F is the functional matrix, and h is the uniform spacing
    
    [gx,gy] = gradient(F, hx, hy);
        
end

function CoF = computeCoF(Fx,Fy, ang)
    
    %First compute the magnitude and direction of contact force
    
    F_contact = [Fx, Fy, 0];
    
    %To ensure contact normal always points outward
    if ang<=0
        contact_normal = [-sin(ang), cos(ang), 0];
    else
        contact_normal = [sin(ang), -cos(ang), 0];
    end

    cone_theta = atan2(norm(cross(F_contact,contact_normal)),dot(F_contact,contact_normal));
    
    CoF = tan(cone_theta);



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


function [U_flex_CU, U_flex_CD] = compute_CUP_CAP_energies(var_theta, var_s)

    smoothed_var_theta = smooth(var_theta, 'rlowess');
    max_logical_vec = islocalmax(smoothed_var_theta).';
    max_indices = find(max_logical_vec == 1, 1, 'first');
    
    if numel(max_indices) == 1
        
        U_flex_CU = computeFlexuralEnergy(...
                            var_theta(1:max_indices),...
                            var_s(1:max_indices));
                        
        U_flex_CD = computeFlexuralEnergy(...
                            var_theta(1+max_indices:end),...
                            var_s(1+max_indices:end));
        
        
    else
        disp('Inflection Issues')
        U_flex_CU = computeFlexuralEnergy(var_theta, var_s);
        U_flex_CD = 0;
    end
    
end



function [X, Y] = sampleEnergyDomain(x_bound, y_bound, size)

    x = linspace(x_bound(1), x_bound(2), size);
    y = linspace(y_bound(1), y_bound(2), size);
    
    [X,Y] = meshgrid(x,y);
end




function stop = outfun(var_theta, optimValues, state)
    stop = false;

    global history;
    if isequal(state,'iter')
          history = [history; var_theta];
    end
        
%     visualize_plots(var_theta)
    
end

function visualize_plots(var_theta)

    global var_s
    
    plot(var_s, var_theta)
    drawnow
    hold on
end
