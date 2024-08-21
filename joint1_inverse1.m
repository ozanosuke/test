function joint1_inverse
clear all
global l1 l2 s1 s2 I1 I2 m1 m2 g v cx cy
l1=1;
l2=1;
s1=0.5;
s2=0.5;
m1=1;
m2=1;
g=9.81;
I1=m1*l1^2/12;
I2=m2*l2^2/12;

v=0.5;
cx=l1*cos(pi/2)+l2*cos(2*pi/3);
cy=l1*sin(pi/2)+l2*sin(2*pi/3);

TSPAN=[0:1E-2:2];
X0=[s1*cos(pi/2) s1*sin(pi/2) pi/2 l1*cos(pi/2)+s2*cos(2*pi/3) l1*sin(pi/2)+s2*sin(2*pi/3) 2*pi/3 0 0 0 0 0 0]';
[t,X]=ode45(@func,TSPAN,X0);

for i=1:length(X)
    figure(1);
    plot([0,0],[0,0]);
    grid;
    axis square;
    xmin=-1.4;xmax=1.4;ymin=-0.4;ymax=2.4;
    axis([xmin,xmax,ymin,ymax]);
    hold on;
    xQ=X(i,1)-s1*cos(X(i,3));
    yQ=X(i,2)-s1*sin(X(i,3));
    xP1=X(i,1)+(l1-s1)*cos(X(i,3));
    yP1=X(i,2)+(l1-s1)*sin(X(i,3));
    xP2=X(i,4)-s2*cos(X(i,6));
    yP2=X(i,5)-s2*sin(X(i,6));
    xR=X(i,4)+(l2-s2)*cos(X(i,6));
    yR=X(i,5)+(l2-s2)*sin(X(i,6));
    plot(xP2,yP2,'o');
    plot([xP1,xQ],[yP1,yQ],'k','linewidth',4);
    plot([xP2,xR],[yP2,yR],'k','linewidth',4);
    drawnow
    hold off
    set(1,'doublebuffer','on')
end

figure(2);
plot(t,X(:,3))
xlabel('time','FontSize',18);
ylabel('X(6)','FontSize',18);


figure(3);
plot(t,X(:,6))
xlabel('time','FontSize',18);
ylabel('phi2','FontSize',18);

function dXdt=func(t,X)
global l1 l2 s1 s2 I1 I2 m1 m2 g v cx cy

M1=[m1 0 0;
    0 m1 0;
    0 0 I1];
M2=[m2 0 0;
    0 m2 0;
    0 0 I2];
Q=[0;
    -m1*g;
    0;
    0;
    -m2*g;
    0];
Cq1=[1 0 s1*sin(X(6));
     0 1 -s1*cos(X(6));
     1 0 -(l1-s1)*sin(X(6));
     0 1  (l1-s1)*cos(X(6));
     0 0 0;
     0 0 0];
  
Cq2=[0 0 0;
    0 0 0;
    -1 0 -s2*sin(X(3));
    0 -1 s2*cos(X(3));
    1 0 -(l2-s2)*sin(X(3));
    0 1 (l2-s2)*cos(X(3))];
A=[M1 zeros(3,3) Cq1';
   zeros(3,3) M2 Cq2';
   Cq1 Cq2 zeros(6,6)];

alph=50;
beta=50;
 
C=[X(4)-s1*cos(X(6));
    X(5)-s1*sin(X(6));
    X(4)+(l1-s1)*cos(X(6))-X(1)+s2*cos(X(3));
    X(5)+(l1-s1)*sin(X(6))-X(2)+s2*sin(X(3));
    X(1)+(l2-s2)*cos(X(3))-(v*t+cx);
    X(2)+(l2-s2)*sin(X(3))-cy];

C1=[X(10)+s1*X(12)*sin(X(6));
    X(11)-s1*X(12)*cos(X(6));
    X(10)-(l1-s1)*X(12)*sin(X(6))-X(7)-s2*X(9)*sin(X(3));
    X(11)+(l1-s1)*X(12)*cos(X(6))-X(8)+s2*X(9)*cos(X(3))
    X(7)-(l2-s2)*sin(X(3))*X(9)-v;
    X(8)+(l2-s2)*cos(X(3))*X(9)];

Gm=-[-s1*X(12)^2*cos(X(6));
    -s1*X(12)^2*sin(X(6));
    (l1-s1)*X(12)^2*cos(X(6))+s2*X(9)^2*cos(X(3));
    (l1-s1)*X(12)^2*sin(X(6))+s2*X(9)^2*sin(X(3));
    -(l2-s2)*cos(X(3))*X(9)^2;
    -(l2-s2)*sin(X(3))*X(9)^2];%-2*alph*C1-(beta^2)*C;
RHS=[Q;
    Gm];
ACC=A\RHS;
dXdt=zeros(12,1);
dXdt(1:6)=X(7:12);
dXdt(7:12)=ACC(1:6);