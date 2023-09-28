function [u_out, theta_out, P_out,...
   theta_id_out, P_id_out,...
   u_buffer_out, y_buffer_out,...
   innov_id_buffer_out, Phibar_out, ubar_out] = DDRCAC_Step(kk, u_in, y_in, z_in, Phi_in,...
                                    P0_val, Rz, Ru, Nf,...
                                    P0_id, N_id, tau_1_id, tau_2_id, eta_id,...
                                    P, theta, Phibar, ubar,...
                                    P_id, theta_id,...
                                    u_buffer, y_buffer, innov_id_buffer)
%% Constants
nu = 1;
Rblock = blkdiag(Rz, Ru);

%% Normalize z
z_in = z_in/(nu+abs(z_in));

%% Update
u_buffer_out = u_buffer;
u_buffer_out(:,2:end) = u_buffer(:,1:end-1);
u_buffer_out(:,1) = u_in;

y_buffer_out = y_buffer;
y_buffer_out(:,2:end) = y_buffer(:,1:end-1);
y_buffer_out(:,1) = y_in;

innov_id_buffer_out = innov_id_buffer;
innov_id_buffer_out(:,2:end) = innov_id_buffer(:,1:end-1);

Phibar_out = [Phi_in; Phibar(1:end-1,:)];

%% DDRCAC
if (kk>1)   
    %% ID stuffs
    ztemp = y_buffer_out(:,2:N_id+1);
    utemp = u_buffer_out(:,1:N_id+1);
    phi_id = kron([-ztemp(:) ; utemp(:)]', eye(1));

    % Update innovation
    innov_id = y_in - phi_id*theta_id;
    innov_id_buffer_out(:,1) = innov_id;

    numTerm = sqrt((1/tau_1_id)*trace(innov_id_buffer(:,1:tau_1_id)*innov_id_buffer(:,1:tau_1_id)'));
    denTerm = sqrt((1/tau_2_id)*trace(innov_id_buffer*innov_id_buffer'));
    if denTerm == 0
        termID = 0;
    else
        termID = numTerm/denTerm ;
    end
    if termID < 1
        termID = 0;
    else
        termID = termID-1;
    end
    betaID = 1 +  eta_id*termID;
    P_id = betaID*P_id;
    P_id =  P_id - P_id*phi_id'/(1 + phi_id*P_id*phi_id')*phi_id*P_id;
    theta_id = theta_id + P_id*phi_id'*(y_in - phi_id*theta_id);

    % Construct Gf
    indexvec = N_id + 1;
    filt_Nf = reshape(-theta_id(indexvec:end), 1, (N_id+1));
    if norm(filt_Nf) == 0 
        filt_Nf(:,1:1) = -ones(1,1);
        disp("Nf set to -1!")
    end

    %% RCAC stuffs
    ubar = [u_in; ubar(1:end-1)];
    u_filt = filt_Nf*ubar;
    Phi_filt = [zeros(1,1) filt_Nf]*Phibar_out;

    % Build Phiblock
    X = [Phi_filt; Phi_in];

    % Compute theta
    gamma = inv(Rblock) + X*P*X';
    P = P - P*X'/gamma*X*P;
    Y = [(Phi_filt*theta + z_in - u_filt); (Phi_in*theta)];
    theta = theta - P*X'*Rblock*Y;
    % theta = theta - P*Phi_filt'*(z_filt - u_filt + Phi_filt*theta) - P*Phi_in'*Rblock(2,2)*Phi_in*theta;


    % ubar = [u_in; ubar(1)];
    % Phibar = [Phi_in; Phibar(1:2,:)];
    % z_filt = z_in;
    % temp = filtNu * ubar;
    % u_filt = temp(1,1);
    % Phi_filt = filtNu * Phibar(2:3,:);
    % 
    % % Build Phiblock
    % X = [Phi_filt; Phi_k];
    % 
    % % Compute theta
    % gamma = inv(Rblock) + X*P*X';
    % P = P - P*X'/gamma*X*P;
    % theta = theta - P*Phi_filt'*(z_filt - u_filt + Phi_filt*theta) - P*Phi_k'*Rblock(2,2)*Phi_k*theta;

else
    P = eye(size(P,1))*P0_val;
    P_id = P0_id;
end

%% Output
temp = Phi_in * theta;
u_out = temp(1,1);

P_out = P;
theta_out = theta;
ubar_out = ubar;

P_id_out = P_id;
theta_id_out = theta_id;



