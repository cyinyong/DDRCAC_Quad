function[linA, linB]=func_linearized_A_B(Jx,Jy,Jz)

%linearized A B
phi=0;the=0;psi=0;
P=0;Q=0;R=0;

A(1,1)=cos(phi)*tan(the)*Q-sin(phi)*tan(the)*R;
A(1,2)=sin(phi)*((sec(the))^2)*Q+cos(phi)*((sec(the))^2)*R;
A(1,4)=1;
A(1,5)=sin(phi)*tan(the);
A(1,6)=cos(phi)*tan(the);

A(2,1)=-sin(phi)*Q-cos(phi)*R;
A(2,5)=cos(phi);
A(2,6)=-sin(phi);

A(3,1)=cos(phi)*sec(the)*Q-sin(phi)*sec(the)*R;
A(3,2)=sin(phi)*sec(the)*tan(the)*Q-cos(phi)*sec(the)*tan(the)*R;
A(3,5)=sin(phi)*sec(the);
A(3,6)=cos(phi)*sec(the);

A(4,5)=((Jy-Jz)*R)/Jx;
A(4,6)=((Jy-Jz)*Q)/Jx;

A(5,4)=((Jz-Jx)*R)/Jy;
A(5,6)=((Jz-Jx)*P)/Jy;

A(6,4)=((Jx-Jy)*Q)/Jz;
A(6,5)=((Jx-Jy)*P)/Jz;

B(4,1)=1/Jx;
B(5,2)=1/Jy;
B(6,3)=1/Jz;

linA=A; 
linB=B;