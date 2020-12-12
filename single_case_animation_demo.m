%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is to present the geometric relation between finger and desk%
% Created by Zachary Chunli JIANG                                         %
% Date: 30.August.2018                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define parameters:
x=-100;%0:10:-100
z=120; %90:1:110
theta=7*pi/180;%0:pi/180:12*pi/180/

p=0.00:0.002:0.3; % pressure is in MPa
gamma=45*pi/180;
m=45;
n=100;
j=0:0.1:18;
k=sqrt(18*18-(j-18).^2);
l=77;
sj=size(j,2);

%point O1:
figure;
x1=x;
z1=z;
pause(1);

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
  
scatter(x1,z1);
axis([-100 100 -100 100],'square');
hold on;
scatter(x2,z2);
hold on;
x0=(0:200)-100;
y0=zeros(1,201);
plot(x0,y0,'k');

for i1=1:size(p,2)
      

      scatter(x5(i1),z5(i1),2,'r');
      %Tip curve
      x6=x4(i1)+j.*sin(phi(i1))-k.*cos(phi(i1));
      z6=z4(i1)-j.*cos(phi(i1))-k.*sin(phi(i1));
      
      Epsilon=0:0.01:delta(i1);
      sE=size(Epsilon,2);
      %Inner surface curve
      x7=x3(i1)+r(i1).*cos(Epsilon+gamma-pi/2-theta);
      z7=z3(i1)-r(i1).*sin(Epsilon+gamma-pi/2-theta);
      
      %External surface curve
      x8=x3(i1)+(r(i1)+22).*cos(Epsilon+gamma-pi/2-theta);
      z8=z3(i1)-(r(i1)+22).*sin(Epsilon+gamma-pi/2-theta);
      
      
      %draw the circle
      th = 0:pi/50:2*pi;
      xcircle=r(i1)*cos(th)+x3(i1);
      ycircle=r(i1)*sin(th)+z3(i1);
      

      line1=plot(x6,z6,'b');
      line2=plot(x7,z7,'b');
      line3=plot(x8,z8,'b');
      line4=plot([x6(sj) x8(sE)],[z6(sj) z8(sE)],'b');
      circle=plot(xcircle,ycircle,'g');
      
      hold on ;
      pause(0.01);
      if i1> 1 && i1<size(p,2)
        %set(line1,'Visible','off');
        set(line2,'Visible','off');
        set(line3,'Visible','off');
        set(line4,'Visible','off');
        set(circle,'Visible','off');
        hold on;
      end
      hold on;
end
  

