% This script is to plot the trajectory of fingertip with optimal grasping
% parameter set
figure;

xo=0;
yo=0;
x0=(0:160)-80;
y0=zeros(1,161);
plot(x0,y0,'k');
hold on;
scatter(xo,yo,'k','filled');
hold on;
axis([-80 80 -20 60]);
hold on;
    
for theta1=0:1:12
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
        %legend('-DynamicLegend');
        scatter(x5(i1),z5(i1),3,'r');
    end
    
    hold all;

end
grid on;
hold off;
