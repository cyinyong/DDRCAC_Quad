%% Compute gains based on poles and damping
clc; clear all; close all
format("default")
addpath Traj_Gen
load gains

%% Set alpha conditions
% alphas = [1.0 1.0].*[0.6:0.2:1.0 1.5:0.5:3.5]';
%alphas = [1.0 1.0].*[0.5:0.1:2.0]';
alphas = [1.0 0.5];

%% Other params
kfwd=0.0;
tau_m = 1/16;
Ka_ol = 1.0;
Ka_il = 0.0;

%% RCAC hyper parameters
rcac_pos_ppi = struct('Ru_xy', 0.005,...
                'Nf_xy', -0.005,...
                'P0_xy', 0.01,...
                'Ru_z', 0.005,...
                'Nf_z', -0.005,...
                'P0_z', 0.005,...
                'PID', 1,...
                'P0_id',0.0005,...
                'N_id', 1,...
                'eta', 0.1,...
                'tau1', 40,...
                'tau2', 200);

rcac_vel_ppi = struct('Ru_xy', 1.0,...
                'Nf_xy', -1.0,...
                'P0_xy', 0.01,...
                'Ru_z', 1.0,...
                'Nf_z', -1.0,...
                'P0_z', 0.001,...
                'PID', 2,...
                'P0_id',0.0005,...
                'N_id', 1,...
                'eta', 0.1,...
                'tau1', 40,...
                'tau2', 200);

rcac_att_ppi = struct('Ru_xy', 1.0,...
                'Nf_xy', -1.0,...
                'P0_xy', 0.1,...
                'Ru_z', 1.0,...
                'Nf_z', -1.0,...
                'P0_z', 0.1,....
                'PID', 1,...
                'P0_id',0.0005,...
                'N_id', 1,...
                'eta', 0.1,...
                'tau1', 40,...
                'tau2', 200);

rcac_rate_ppi = struct('Ru_xy', 1.0,...
                'Nf_xy', -1.0,...
                'P0_xy', 0.000001,...
                'Ru_z', 1.0,...
                'Nf_z', -1.0,...
                'P0_z', 0.000001,...
                'PID', 2,...
                'P0_id',0.0005,...
                'N_id', 1,...
                'eta', 0.1,...
                'tau1', 40,...
                'tau2', 200);

%% Trajectory generation
randn('state',1)
rand('state',1)

a_max   = 1;
v_max   = 2;

Waypoints = readmatrix('Trajectory_Hilbert.xlsx');

wait_time = 2;      % wait at a waypoint before starting
[Traj,time] = func_Stitch_trajectory(Waypoints,a_max, v_max,wait_time,time_step);

simin = timeseries(Traj,time);
T = numel(time)*time_step;

%% Create bus object
clear elems
kk = 1;
elems(kk) = Simulink.BusElement;
elems(kk).Name = 'Ru_xy'; kk=kk+1;
elems(kk) = Simulink.BusElement;
elems(kk).Name = 'Nf_xy'; kk=kk+1;
elems(kk) = Simulink.BusElement;
elems(kk).Name = 'P0_xy'; kk=kk+1;
elems(kk) = Simulink.BusElement;
elems(kk).Name = 'Ru_z'; kk=kk+1;
elems(kk) = Simulink.BusElement;
elems(kk).Name = 'Nf_z'; kk=kk+1;
elems(kk) = Simulink.BusElement;
elems(kk).Name = 'P0_z'; kk=kk+1;
elems(kk) = Simulink.BusElement;
elems(kk).Name = 'PID'; kk=kk+1;
elems(kk) = Simulink.BusElement;
elems(kk).Name = 'P0_id'; kk=kk+1;
elems(kk) = Simulink.BusElement;
elems(kk).Name = 'N_id'; kk=kk+1;
elems(kk) = Simulink.BusElement;
elems(kk).Name = 'eta'; kk=kk+1;
elems(kk) = Simulink.BusElement;
elems(kk).Name = 'tau1'; kk=kk+1;
elems(kk) = Simulink.BusElement;
elems(kk).Name = 'tau2'; kk=kk+1;
RCACParamsBus = Simulink.Bus;
RCACParamsBus.Elements = elems;

%% Do simulations

for i = 1:length(alphas(:,1))
    alpha_inner = alphas(i,1);
    alpha_outer = alphas(i,2);

    % out_fsb(i) = sim('FSF_pole_placement_loops_V5.slx',T);
    out_ppi(i) = sim('PPI_loops.slx',T);

    if out_ppi(i).tout(end) == T
        pos_mse_ppi(i) = mean(vecnorm(out_ppi(i).err_pos,2,2));
        eul_mse_ppi(i) = mean(vecnorm(out_ppi(i).err_eul,2,2));
    else
        pos_mse_ppi(i) = 3.0;
        eul_mse_ppi(i) = 3.0;
    end
end

%% Plot
% Which alpha
ii = 1;

figure;
plot3(out_ppi(ii).ref_pos(:,1),out_ppi(ii).ref_pos(:,2),out_ppi(ii).ref_pos(:,3),...
    out_ppi(ii).pos(:,1),out_ppi(ii).pos(:,2),out_ppi(ii).pos(:,3))
grid
title('PPI Traj')


figure
subplot(313)
plot(out_ppi(ii).tout,out_ppi(ii).ref_eul(:,3),out_ppi(ii).tout,out_ppi(ii).eul(:,3))
grid
title('PPI Heading')

subplot(311)
plot(out_ppi(ii).tout,out_ppi(ii).ref_eul(:,1),out_ppi(ii).tout,out_ppi(ii).eul(:,1))
grid
title('PPI Phi')

subplot(312)
plot(out_ppi(ii).tout,out_ppi(ii).ref_eul(:,2),out_ppi(ii).tout,out_ppi(ii).eul(:,2))
grid
title('PPI Theta')


%% 
figure
subplot(311)
plot(out_ppi(ii).tout,out_ppi(ii).rcac_vel_u(:,1),out_ppi(ii).tout,out_ppi(ii).vel_u(:,1))
grid
subplot(312)
plot(out_ppi(ii).tout,out_ppi(ii).rcac_vel_u(:,2),out_ppi(ii).tout,out_ppi(ii).vel_u(:,2))
grid
subplot(313)
plot(out_ppi(ii).tout,out_ppi(ii).rcac_vel_u(:,3),out_ppi(ii).tout,out_ppi(ii).vel_u(:,3))
grid
