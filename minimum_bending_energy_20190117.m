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


figure;

%%
global options;
options = optimoptions('fmincon', 'OutputFcn', @outfun,'Display','iter-detailed','Algorithm','interior-point');
options.MaxFunctionEvaluations = 1000000;
options.MaxIterations = 3000;
options.OptimalityTolerance=1e-3;
options.ConstraintTolerance=1e-3;
options.StepTolerance=1e-8;

%%
% Invoke optimization with a strating guess
global steps;
global s L;
global x_end y_end;
global concave_up_flex concave_down_flex total_flex_energy;
global collect_x_end collect_y_end;
global collect_a;
%s is in terms of the coordinate system of the object
%L is the total length of the paper strip
steps=100;
L=250; 
s=linspace(0,L,steps);

%%
xo=0;
yo=0;
x0=(0:(L*1.2+160))-L*1.2-80;
y0=zeros(1,(L*1.2+161));
plot(x0,y0,'k');
axis([-L*1.2-80 80 -20 200]);
drawnow
hold on;
drawnow
scatter(xo,yo,'k','filled');
drawnow
hold on;



%%
%This is the main logic
scale=0.05;%% minimal 0.02

for ii=0:scale:0.5
    for jj=0:scale:0.5
        %if jj<ii
            x_end=L-ii*L;           
            y_end=jj*L;
            if sqrt(x_end^2+y_end^2)<=L
                start_optimization();
            end
        %end
    end
end
grid on;
hold off;
%%
% plot out the flex energy
% figure;
% nof=1:1:size(concave_down_flex,2);
% plot(nof,concave_up_flex);
% hold on;
% plot(nof,concave_down_flex);
% hold on;
% plot(nof,total_flex_energy);
% legend('concave up', 'concave down', 'total flex energy');
% hold off;

smap=size(concave_down_flex,2);
energy_map=zeros(smap,6);

figure
for sx=1:smap
    energy_map(sx,1)=collect_x_end(sx);
    energy_map(sx,2)=collect_y_end(sx);
    energy_map(sx,3)=concave_up_flex(sx);
    energy_map(sx,4)=concave_down_flex(sx);
    energy_map(sx,5)=total_flex_energy(sx);
    
    current_a=collect_a((sx*10-9):(sx*10));
    k=curvature(current_a);
    energy_map(sx,6)=max(k);
    plot(s(1:steps-1)/L,k);
    hold on;
end
hold off;

figure;
[qx,qy]=meshgrid(L*0.5:L*scale:L*1,0:L*scale:L*0.5);
F1=TriScatteredInterp(energy_map(:,1),energy_map(:,2),energy_map(:,3));
qz=F1(qx,qy);
surf(qx,qy,qz);
view(2);
colorbar;
hold off;

[px,py] = gradient(qz, scale, scale); 

starty = 0:L*scale:L*0.5;
startx = 0.5*L*ones(size(starty));
figure 
contour(qx,qy,qz);
hold on;
quiver(qx,qy,-px,-py);
hold on;
streamline(qx,qy,-px,-py,startx,starty);
hold off;


figure
F2=TriScatteredInterp(energy_map(:,1),energy_map(:,2),energy_map(:,6));
cz=F2(qx,qy);
surf(qx,qy,cz);
view(2);
colorbar;
hold off;

[cx,cy] = gradient(cz, scale, scale);
figure 
contour(qx,qy,cz);
hold on;
quiver(qx,qy,cx,cy);
hold on;
streamline(qx,qy,cx,cy,startx,starty);
hold off;
%%
function start_optimization()

    global concave_up_flex concave_down_flex total_flex_energy;
    global collect_x_end collect_y_end;
    global collect_a;
    global s L options;
    global x_end y_end;
    a=[0.2,0.2,0.6,0.6,0.1,0.1,0.1,0.1,0.1,0.1];
    [a,fval] = fmincon(@total_energy, a,[],[],[],[],[],[], @constrain_function, options);
    theta=ritz_fourier_model(a);

    xc = L*cumtrapz(s/L, cos(theta))-L;
    yc = L*cumtrapz(s/L, sin(theta));

    plot(xc, yc);
    drawnow
    hold on;
    
    %calculate flex energies;

    if concave_up_flex_energy(a)>=0 && concave_down_flex_energy(a)>=0
        concave_up_flex=[concave_up_flex, concave_up_flex_energy(a)];
        concave_down_flex=[concave_down_flex, concave_down_flex_energy(a)];
        total_flex_energy=concave_up_flex+concave_down_flex;
        collect_x_end=[collect_x_end,x_end];
        collect_y_end=[collect_y_end,y_end];
        collect_a=[collect_a,a]
    end

%     disp('concave_up_energy:');
%     disp(concave_up_flex);
%     disp('concave_down_energy:');
%     disp(concave_down_flex);
%     disp('total_flex_energy:');
%     disp(total_flex_energy);
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
