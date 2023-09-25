clear all; close all;

%% Properties
% Physical properties
m = 0.8;
Jx = 0.005;Jy = 0.005;Jz = 0.009;
g = [0; 0; 10];

time_step = 0.001;   % time step of the trajectory

%% Controller design targets 
% w_vxy = 1.0;
% w_xy = 0.2;
% w_vz = 2.5;
% w_z = 0.5; 
% w_pq = 8.0;
% w_tet = 2.0;
% w_r = 3.0;
% w_psi = 0.5;

w_vxy = 1.0;
w_xy = 0.2;
w_vz = 2.5;
w_z = 0.5; 
w_pq = 8.0;
w_tet = 2.0;
w_r = 3.0;
w_psi = 0.5;
zeta_ppi = 0.8;

% p_pq = [-11.0 -11.1 -11.2];
p_pq = [-11.0 -11.1 -11.2];
p_r = [-3.0 -3.1 -3.2];
p_xy = [-1.8 -1.7 -1.6];
p_z = [-2.0 -2.1 -2.2];

%% Outer loop dynamics
A = [zeros(3,3),eye(3);zeros(3,3),zeros(3,3)];
B = [zeros(3,3);diag([1/m,1/m,1/m])];
C = eye(6);
C_r = C(1:3,:);     % Position output
C_v = C(4:6,:);     % Velocity Output

% Compute linearized models in tf
linC = [zeros(3,3) eye(3)];
linD = zeros(3,3);
[b,a] = ss2tf(A,B,linC,linD,1);
P_lin = minreal(tf(b(1,:),a));
P_lin_d = c2d(P_lin,time_step);

% P-PI Controller for XY axis
[kp1, kp2, ki2, pole_xy] = compute_ppi_gains_discrete(P_lin(1),w_xy,w_vxy,zeta_ppi,time_step);
Kp1_OL = eye(3)*kp1;  
Kp2_OL = eye(3)*kp2;
Ki2_OL = eye(3)*ki2;

% P-PI Controller for Z axis
[kp1, kp2, ki2, pole_z] = compute_ppi_gains_discrete(P_lin(1),w_z,w_vz,zeta_ppi,time_step);
Kp1_OL(3,3) = kp1;  
Kp2_OL(3,3) = kp2;
Ki2_OL(3,3) = ki2;

% FSB controller
A_single = [0 1 0; 0 0 0; -1 0 0];
B_single = [0; B(1+3,1); 0];
sysd = c2d(ss(A_single,B_single,eye(3),0),time_step);
for i=1:3
    if i < 3
        [kfsb(i,:), pole_outer_fsb(i,:)] = compute_fsb_gains_v2(sysd.A,sysd.B,exp(p_xy.*time_step));
    else
        [kfsb(i,:), pole_outer_fsb(i,:)] = compute_fsb_gains_v2(sysd.A,sysd.B,exp(p_z.*time_step));
    end
end

% Pole placement
K_o     = [diag(-kfsb(:,1)) diag(-kfsb(:,2)) diag(-kfsb(:,3))];
Kx_o    = K_o(:,1:6);
Kq_o    = K_o(:,7:9);

%% Inner loop dynamics

C_q = C_r;
C_w = C_v;

% Full state feedback
[linA, linB]    = func_linearized_A_B(Jx,Jy,Jz);

for i=1:3
    linC = [zeros(3,3) eye(3)];
    linD = zeros(3,3);
    [b,a] = ss2tf(linA,linB,linC,linD,i);
    P_lin(i) = minreal(tf(b(i,:),a)); 

    if i<3
        [kr(i), kv(i), ki(i), pole_inner_ppi(i,:)] = compute_ppi_gains_discrete(P_lin(i),w_tet,w_pq,zeta_ppi,time_step);
    else
        [kr(i), kv(i), ki(i), pole_inner_ppi(i,:)] = compute_ppi_gains_discrete(P_lin(i),w_psi,w_r,zeta_ppi,time_step);
    end
end

Kp1_IL = diag(kr);
Kp2_IL = diag(kv);
Ki2_IL = diag(ki);

% Compute FSB controller poles
for i=1:3
    if i<3
        B_single = [0; linB(i+3,i); 0];
        sysd = c2d(ss(A_single,B_single,eye(3),0),time_step);
        [kfsb(i,:), pole_inner_fsb(i,:)] = compute_fsb_gains_v2(sysd.A,sysd.B,exp(p_pq.*time_step));
    else
        [kfsb(i,:), pole_inner_fsb(i,:)] = compute_fsb_gains_v2(sysd.A,sysd.B,exp(p_r.*time_step));
    end
end

K_i     = [diag(-kfsb(:,1)) diag(-kfsb(:,2)) diag(-kfsb(:,3))];
Kx_i    = K_i(:,1:6);
Kq_i    = K_i(:,7:9);

%% Save controller
save("gains.mat","Kx_i","Kq_i","Kx_o","Kq_o","Kp1_OL","Kp2_OL","Ki2_OL","Kp1_IL","Kp2_IL","Ki2_IL","time_step","Jx","Jy","Jz","m","C_q","C_r","C_v","C_w");
