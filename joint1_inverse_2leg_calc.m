syms x1 x2 x3 x4 y1 y2 y3 y4 phi1 phi2 phi3 phi4 l1 l2 l3 l4 s1 s2 s3 s4 cx1 cx2 cy1 cy2 v t w X1 X2 X3 X4 Y1 Y2 Y3 Y4 PHI1 PHI2 PHI3 PHI4
C=[x1-s1*cos(phi1);
    y1-s1*sin(phi1);
    x1+(l1-s1)*cos(phi1)-x2+s2*cos(phi2);
    y1+(l1-s1)*sin(phi1)-y2+s2*sin(phi2);
    x2+(l2-s2)*cos(phi2)-(v*t+cx1);
    y2+(l2-s2)*sin(phi2)-cy1;
    x2+(l2-s2)*cos(phi2)-x4-(l4-s4)*cos(phi4);
    y2+(l2-s2)*sin(phi2)-y4-(l4-s4)*sin(phi4);   
    x3+(l3-s3)*cos(phi3)-x4-s4*cos(phi4);
    y3+(l3-s3)*sin(phi3)-y4-s4*sin(phi4);
    x3+s3*cos(phi3)-(w*t+cx2);
    y3+s3*sin(phi3)-cy2];

Cq=jacobian(C,[x1,y1,phi1,x2,y2,phi2,x3,y3,phi3,x4,y4,phi4])
dq=[X1; Y1; PHI1;X2;Y2;PHI3;X3;Y3;PHI3;X4;Y4;PHI4];
C1=Cq*dq
C1q=jacobian(C1,[x1,y1,phi1,x2,y2,phi2,x3,y3,phi3,x4,y4,phi4])
C2=C1q*dq