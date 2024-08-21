syms l1 l2 m1 m2 I1 I2 x1 y1 phi1 x2 y2 phi2 X1 Y1 PHI1 X2 Y2 PHI2 v t cx cy

C=[x1-l1/2*cos(phi1);
    y1-l1/2*sin(phi1);
    x1+l1/2*cos(phi1)-x2+l2/2*cos(phi2);
    y1+l1/2*sin(phi1)-y2+l2/2*sin(phi2);
    x2+l2/2*cos(phi2)-(v*t+cx);
    y2+l2/2*sin(phi2)-cy];

Cq1=jacobian(C,[x1,y1,phi1])
    
Cq2=jacobian(C,[x2,y2,phi2])

Cq=cat(2,Cq1,Cq2);
q1=[X1;Y1;PHI1;X2;Y2;PHI2];
C1=Cq*q1
C1q=jacobian(C1,[x1,y1,phi1,x2,y2,phi2]);
C2=C1q*q1
Ct=jacobian(C,t);
Ctq=jacobian(Ct,[x1,y1,phi1,x2,y2,phi2])
