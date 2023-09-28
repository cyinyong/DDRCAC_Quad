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
rcac_pos_ppi = struct('Ru_xy', 0.01,...
                'Nf_xy', -0.01/2,...
                'P0_xy', 0.1,...
                'Ru_z', 0.01,...
                'Nf_z', -0.01/2,...
                'P0_z', 0.1,...
                'PID', 1,...
                'P0_id',1,...
                'N_id', 1,...
                'eta', 0.1,...
                'tau1', 15,...
                'tau2', 20);

rcac_vel_ppi = struct('Ru_xy', 0.01,...
                'Nf_xy', -1.0/2,...
                'P0_xy', 1.0,...
                'Ru_z', 0.01,...
                'Nf_z', -1.0/2,...
                'P0_z', 1.0,...
                'PID', 2,...
                'P0_id',1,...
                'N_id', 1,...
                'eta', 0.1,...
                'tau1', 15,...
                'tau2', 20);

rcac_att_ppi = struct('Ru_xy', 1.0,...
                'Nf_xy', -1.0,...
                'P0_xy', 0.1,...
                'Ru_z', 1.0,...
                'Nf_z', -1.0,...
                'P0_z', 0.1,....
                'PID', 1,...
                'P0_id',0.0005,...
                'N_id', 1,...
                'eta', 0.05,...
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

    out_ppi(i) = sim('PPI_loops.slx',T);
    out_rcac(i) = sim('PPI_loops_rcac.slx',T);
    out_ddrcac(i) = sim('PPI_loops_ddrcac.slx',T);
    
    rcac_pos_ppi.Nf_xy = -rcac_pos_ppi.Nf_xy ;
    rcac_vel_ppi.Nf_xy = -rcac_vel_ppi.Nf_xy ;
    % rcac_pos_ppi.Nf_z = 0.1*rcac_pos_ppi.Nf_z ;
    % rcac_vel_ppi.Nf_z = 0.1*rcac_vel_ppi.Nf_z ;
    out_rcac2(i) = sim('PPI_loops_rcac.slx',T);

    if out_ppi(i).tout(end) == T
        pos_mse(i,1) = mean(vecnorm(out_ppi(i).err_pos,2,2));
    else
        pos_mse(i,1) = 3.0;
    end

    if out_ddrcac(i).tout(end) == T
        pos_mse(i,2) = mean(vecnorm(out_ddrcac(i).err_pos,2,2));
    else
        pos_mse(i,2) = 3.0;
    end

    if out_rcac(i).tout(end) == T
        pos_mse(i,3) = mean(vecnorm(out_rcac(i).err_pos,2,2));
    else
        pos_mse(i,3) = 3.0;
    end

    if out_rcac2(i).tout(end) == T
        pos_mse(i,4) = mean(vecnorm(out_rcac2(i).err_pos,2,2));
    else
        pos_mse(i,4) = 3.0;
    end
      
        
end

%% Plot
% Which alpha
ii = 1;


%%
close all;
% h=figure(1);
% set(0,'defaultAxesFontName', 'times',...
%     'defaultTextFontName','times',...
%     'DefaultAxesFontSize', 10, ...
%     'DefaultTextFontSize', 10, ...
%     'DefaultLineLineWidth',1,...
%     'DefaultLineMarkerSize', 1,...
%     'DefaultTextInterpreter','latex')
% h.Color='w';
% h.Units = 'inches';
% % h.Position=[.25 .25 3.85 6.3];
% % h.PaperPosition=[.25 .25 3.85 6.3];
% ha = tight_subplot(1,1,[.02 0.02],[.1 .05],[.1 .05]);
% fg=[0 102 51]/255;
% fr=[1 .2 .2];
% CO= [  1 .2 .2
%      0 .6 .3       
%     0    0.4470    0.7410
% 
%     0.9290    0.6940    0.1250
%     0.4940    0.1840    0.5560
%     0.4660    0.6740    0.1880
%     0.3010    0.7450    0.9330
%     0.6350    0.0780    0.1840];
% set(groot,'defaultAxesColorOrder',CO)
% 
% axes(ha(1))
% ppp = plot3(out_ppi(ii).pos(:,1),out_ppi(ii).pos(:,2),out_ppi(ii).pos(:,3),...
%         out_ddrcac(ii).pos(:,1),out_ddrcac(ii).pos(:,2),out_ddrcac(ii).pos(:,3),...
%         out_rcac(ii).pos(:,1),out_rcac(ii).pos(:,2),out_rcac(ii).pos(:,3),...
%         out_rcac2(ii).pos(:,1),out_rcac2(ii).pos(:,2),out_rcac2(ii).pos(:,3));
% hold on
% p1 = plot3(out_ppi(ii).ref_pos(:,1),out_ppi(ii).ref_pos(:,2),out_ppi(ii).ref_pos(:,3),'k--');
% xlabel('$\bf x(m)$')
% ylabel('$\bf y(m)$')
% zlabel('$\bf z(m)$')
% xlim([-2 8])
% ylim([-8 2])
% 
% h1 = legend([p1],'$Ref$','box','off');
% set(h1,'Interpreter','latex','location','northeast')
% a=axes('position',get(gca,'position'),'visible','off');
% h2 = legend(a,[ppp],'PPI','DDRCAC','RCAC, $G_{\rm f}=-{\beta}/{z}$','RCAC, $G_{\rm f}={\beta}/{z}$');
% set(h2,'Interpreter','latex','location','south')


h=figure(2);
set(0,'defaultAxesFontName', 'times',...
    'defaultTextFontName','times',...
    'DefaultAxesFontSize', 10, ...
    'DefaultTextFontSize', 10, ...
    'DefaultLineLineWidth',1,...
    'DefaultLineMarkerSize', 1,...
    'DefaultTextInterpreter','latex')
h.Color='w';
h.Units = 'inches';

ha = tight_subplot(1,1,[.02 0.02],[.1 .05],[.1 .05]);
fg=[0 102 51]/255;
fr=[1 .2 .2];
CO= [  1 .2 .2
     0 .6 .3       
    0    0.4470    0.7410
    
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];
set(groot,'defaultAxesColorOrder',CO)

axes(ha(1))
ppp = plot(out_ppi(ii).pos(:,1),out_ppi(ii).pos(:,2),...
        out_ddrcac(ii).pos(:,1),out_ddrcac(ii).pos(:,2),...
        out_rcac(ii).pos(:,1),out_rcac(ii).pos(:,2),...
        out_rcac2(ii).pos(:,1),out_rcac2(ii).pos(:,2));
hold on
p1 = plot(out_ppi(ii).ref_pos(:,1),out_ppi(ii).ref_pos(:,2),'k--');
xlabel('$\bf x(m)$')
ylabel('$\bf y(m)$')
xlim([-2 8])
ylim([-8 2])

h1 = legend([p1],'$Ref$','box','off');
set(h1,'Interpreter','latex','location','northwest')
a=axes('position',get(gca,'position'),'visible','off');
h2 = legend(a,[ppp],'PPI','DDRCAC','RCAC, $G_{\rm f}=-{\beta}/{z}$','RCAC, $G_{\rm f}={\beta}/{z}$');
set(h2,'Interpreter','latex','location','southwest')

h=figure(3);
set(0,'defaultAxesFontName', 'times',...
    'defaultTextFontName','times',...
    'DefaultAxesFontSize', 10, ...
    'DefaultTextFontSize', 10, ...
    'DefaultLineLineWidth',1,...
    'DefaultLineMarkerSize', 1,...
    'DefaultTextInterpreter','latex')
h.Color='w';
h.Units = 'inches';
ha = tight_subplot(3,1,[.02 0.02],[.1 .05],[.1 .05]);
fg=[0 102 51]/255;
fr=[1 .2 .2];
CO= [  1 .2 .2
     0 .6 .3       
    0    0.4470    0.7410
    
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];
set(groot,'defaultAxesColorOrder',CO)

for kk =1:3
    axes(ha(kk))
    ppp = plot(out_ppi(ii).tout,out_ppi(ii).pos(:,kk),...
            out_ddrcac(ii).tout,out_ddrcac(ii).pos(:,kk),...
            out_rcac(ii).tout,out_rcac(ii).pos(:,kk),...
            out_rcac2(ii).tout,out_rcac2(ii).pos(:,kk));
    hold on
    p1 = plot(out_ppi(ii).tout,out_ppi(ii).ref_pos(:,kk),'k--');
    xlim([0 70])
    axP = get(gca,'Position');
    set(gca,'Position',axP)
    if kk==1
        ylabel('$\bf x(m)$')
        set(gca,'xticklabel',[])
        ylim([0 6])
    elseif kk==2
        ylabel('$\bf y(m)$')
        set(gca,'xticklabel',[])
        ylim([-6 0])
    else
        ylabel('$\bf z(m)$')
        xlabel('$\bf t(s)$')
        ylim([-2 6])
    end
end
%%
h = figure(4);
set(0,'defaultAxesFontName', 'times',...
    'defaultTextFontName','times',...
    'DefaultAxesFontSize', 10, ...
    'DefaultTextFontSize', 10, ...
    'DefaultLineLineWidth',1,...
    'DefaultLineMarkerSize', 1,...
    'DefaultTextInterpreter','latex')
h.Color='w';
h.Units = 'inches';
ha = tight_subplot(1,1,[.02 0.02],[.22 .05],[.1 .05]);
fg=[0 102 51]/255;
fr=[1 .2 .2];
CO= [  1 .2 .2
     0 .6 .3       
    0    0.4470    0.7410
    
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];
set(groot,'defaultAxesColorOrder',CO)

axes(ha(1))
b=bar(pos_mse(ii,:));
b.FaceColor = 'flat';
for jj = 2:4
    b.CData(jj,:) = CO(jj,:);
end
set(gca,'XTickLabel',{'\bf{PPI}','\bf{DDRCAC}','\bf{RCAC, $G_{\rm f}=-{\beta}/{z}$}','\bf{RCAC, $G_{\rm f}={\beta}/{z}$}'},'TickLabelInterpreter','latex');
ylabel('$Pos M.S.E(m)$')

%% 

for jj = 1:length(out_ddrcac(ii).tout)
    pos_norm(jj,1) = norm(out_ddrcac(ii).theta_id_pos(jj,1:3));
    pos_norm(jj,2) = norm(out_ddrcac(ii).theta_id_pos(jj,4:6));
    pos_norm(jj,3) = norm(out_ddrcac(ii).theta_id_pos(jj,7:9));

    vel_norm(jj,1) = norm(out_ddrcac(ii).theta_id_vel(jj,1:3));
    vel_norm(jj,2) = norm(out_ddrcac(ii).theta_id_vel(jj,4:6));
    vel_norm(jj,3) = norm(out_ddrcac(ii).theta_id_vel(jj,7:9));
end

figure
subplot(3,1,1)
plot(out_ddrcac(ii).tout,out_ddrcac(ii).theta_id_pos(:,1:3),out_ddrcac(ii).tout,pos_norm(:,1),'k--')
ylabel('x axis')
subplot(3,1,2)
plot(out_ddrcac(ii).tout,out_ddrcac(ii).theta_id_pos(:,4:6),out_ddrcac(ii).tout,pos_norm(:,2),'k--')
ylabel('y axis')
subplot(3,1,3)
plot(out_ddrcac(ii).tout,out_ddrcac(ii).theta_id_pos(:,7:9),out_ddrcac(ii).tout,pos_norm(:,3),'k--')
ylabel('z axis')
legend('thetaID[0]','thetaID[1]','thetaID[2]','norm(theta)','location','best')
sgtitle('ThetaID for Position P controller loop')

figure
subplot(3,1,1)
plot(out_ddrcac(ii).tout,out_ddrcac(ii).theta_id_vel(:,1:3),out_ddrcac(ii).tout,vel_norm(:,1),'k--')
ylabel('x axis')
subplot(3,1,2)
plot(out_ddrcac(ii).tout,out_ddrcac(ii).theta_id_vel(:,4:6),out_ddrcac(ii).tout,vel_norm(:,2),'k--')
ylabel('y axis')
subplot(3,1,3)
plot(out_ddrcac(ii).tout,out_ddrcac(ii).theta_id_vel(:,7:9),out_ddrcac(ii).tout,vel_norm(:,3),'k--')
ylabel('z axis')
legend('thetaID[0]','thetaID[1]','thetaID[2]','norm(theta)','location','best')
sgtitle('ThetaID for velocity PI controller loop')