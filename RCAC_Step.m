function [u_out, theta_out, P_out, ubar_out, Phibar_out] = RCAC_Step(kk, u_in, z_in, phi_in, P0_val, Rz, Ru, Nf, errorNormMode,...
                                        P, theta, ubar, Phibar)
%% Constants
mu = 1.0;
nu = 1.0;

% nf = 2;
% lu = 1;
% lz = 1;

filtNu = [0 Nf];
Rblock = blkdiag(Rz, Ru);

%% Set inputs for current step
Phi_k = phi_in;
%% Normalize z
z_k = z_in/(nu+abs(z_in));
u_k = u_in;

%% Normalize error
if errorNormMode > 1
    switch errorNormMode
        case 1
            z_k = mu*nu*z_k / (mu + nu *abs(z_k));
        case 2
            z_k = (2*nu/pi)*atan(pi*nu*z_k)/2/mu;
        case 3
            z_k = (mu*nu*z_k)/sqrt(mu*mu + nu*nu*z_k*z_k);
        case 4
            z_k = mu*tanh((nu*z_k)/mu);
        case 5
            z_k = mu*erf((sqrt(pi)*nu*z_k)/(2*mu));
    end
end

%% Compute theta
if (kk>3)
    % Filter data
    ubar = [u_k; ubar(1)];
    Phibar = [Phi_k; Phibar(1:2,:)];
    z_filt = z_k;
    temp = filtNu * ubar;
    u_filt = temp(1,1);
    Phi_filt = filtNu * Phibar(2:3,:);

    % Build Phiblock
    X = [Phi_filt; Phi_k];
 
    % Compute theta
    gamma = inv(Rblock) + X*P*X';
    P = P - P*X'/gamma*X*P;
    theta_out = theta - P*Phi_filt'*(z_filt - u_filt + Phi_filt*theta) - P*Phi_k'*Rblock(2,2)*Phi_k*theta;
else
    theta_out = theta;
    P = eye(size(P,1))*P0_val;
end

%% Compute outputs
temp = Phi_k * theta;
u_out = temp(1,1);

P_out = P;
ubar_out = ubar;
Phibar_out = Phibar;

end