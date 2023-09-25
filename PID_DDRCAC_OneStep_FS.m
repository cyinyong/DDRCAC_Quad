function [u_out, theta_out, P_out, theta_id_out, P_id_out, lambda_id_out] = PID_DDRCAC_OneStep_FS(kk, u_in, y_in, z_in, kdel, N_id, FLAG, FREEZE)
persistent ltheta ltheta_id lphi lu ly lz Rbar %Constants
persistent theta_buffer theta_id_buffer u_buffer y_buffer z_buffer integral_z internal_u_buffer innovID_buffer %Data buffers
persistent PHI PHI_bar u_bar PHI_filt u_filt  % Regressors
persistent P lambda_inv theta_0 % Recursive Variables
persistent P_id % Recursive Variables for Identification
persistent tau_1_id tau_2_id eta_id

if kk == 1
    lu      = size(u_in(:,1),1);
    ly      = size(y_in,1);
    lz      = size(z_in,1);
    ltheta  = 0;
    lphi    = zeros(lu,lz);
    for ii = 1:lu
        for jj = 1:lz
            if FLAG.PID(ii,jj) == 1
                ltheta = ltheta + 1;
                lphi(ii,jj)   = 1;
            elseif FLAG.PID(ii,jj) == 2
                ltheta = ltheta + 2;
                lphi(ii,jj)   = 2;
            elseif FLAG.PID(ii,jj) == 3
                ltheta = ltheta + 3;
                lphi(ii,jj)   = 3;
            end
        end
    end
    Nc_max      = max(abs(FLAG.PID(:)));
    
    ltheta_id = N_id*ly*ly + (N_id+1)*ly*lu;
    theta_id_buffer = zeros(ltheta_id,2);

    lambda_inv  = 1/FLAG.lambda;
    
    theta_0             = zeros(ltheta,1);
    theta_buffer        = [theta_0 theta_0];
    u_buffer            = zeros(lu,max(Nc_max,N_id)+2); %[ u(k-1) ...]
    internal_u_buffer   = zeros(lu,lz,Nc_max+1); %[ u_i(k-1) ...]
    y_buffer            = zeros(ly,max(Nc_max,N_id)+2); %[ z(k-1) ...]
    z_buffer            = zeros(lz,max(Nc_max,N_id)+2); %[ z(k-1) ...]
    innovID_buffer      = zeros(ly,FLAG.tau_2);

    integral_z          = zeros(lz,1);
            
    PHI = zeros(lu,ltheta); %Phi
    PHI_bar = zeros((N_id+2)*lu,ltheta);
    PHI_filt = zeros(lz,ltheta);
    u_bar   = zeros((N_id+1)*lu,1);
    u_filt   = zeros(lz,1);
    P        = FLAG.P0;
    Rbar = blkdiag(FLAG.Rz,FLAG.Ru);

    P_id = FLAG.P0_id;
    tau_1_id = FLAG.tau_1;
    tau_2_id = FLAG.tau_2;
    eta_id = FLAG.eta;
   
    u_out = zeros(lu,1);
    theta_out = theta_buffer(:,1);
    P_out = P;
    theta_id_out = theta_id_buffer(:,1);
    P_id_out = P_id;
    lambda_id_out = 0;

else
    u_buffer(:,2:end) = u_buffer(:,1:end-1);
    y_buffer(:,2:end) = y_buffer(:,1:end-1);
    z_buffer(:,2:end) = z_buffer(:,1:end-1);
    internal_u_buffer(:,:,2:end) = internal_u_buffer(:,:,1:end-1);
    innovID_buffer(:,2:end) = innovID_buffer(:,1:end-1);
    if kk>kdel
        u_buffer(:,1) = u_in;   %u(k-1)
        y_buffer(:,1) = y_in;   %y(k-1)
        z_buffer(:,1) = z_in;   %z(k-1)
        %No need to add data to internal_u_buffer yet, will be done after
        %constructing PHI
    end
    
    %%%%%%%%% ID Stuff %%%%%%%%%%%%%

    ztemp = y_buffer(:,2:N_id+1);
    utemp = u_buffer(:,1:N_id+1);
    phi_ID = kron( [ -ztemp(:) ; utemp(:) ]' , eye(ly) );

    innovID = y_in - phi_ID*theta_id_buffer(:,1);
    innovID_buffer(:,1) = innovID;

    numTerm = sqrt( (1/tau_1_id)*trace ( innovID_buffer(:,1:tau_1_id)*innovID_buffer(:,1:tau_1_id)' ) );
    denTerm = sqrt( (1/tau_2_id)*trace ( innovID_buffer*innovID_buffer' ) );
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
    

    if FREEZE
        ThetaID = theta_id_buffer(:,1);
    else
        P_id = betaID*P_id;
        P_id =  P_id - P_id*phi_ID'/( eye(ly) + phi_ID*P_id*phi_ID')*phi_ID*P_id;
        ThetaID = theta_id_buffer(:,1) + P_id*phi_ID'*( y_in - phi_ID*theta_id_buffer(:,1) );
        theta_id_buffer(:,2) = theta_id_buffer(:,1);
        theta_id_buffer(:,1) = ThetaID;
    end
    
    %%% Construct Gf %%%
    indexvec = (N_id*ly*ly) + 1;
    filt_Nf = reshape( -theta_id_buffer(indexvec:end,1) , ly , (N_id+1)*lu );
    if norm(filt_Nf) == 0
        filt_Nf(:,1:lu) = -ones(ly,lu);
    end

    %%%%%%%%% RCAC Stuff %%%%%%%%%%%%%
    
    integral_z = integral_z + z_buffer(:,2);
    %Construct PHI
    count = 1;
    for iii = 1:lu
        count_internal = 1;
        phi_temp = zeros(1,sum(lphi(iii,:)));
        for jjj = 1:lz
                if FLAG.PID(iii,jjj) == 1
                    phi_temp(:,count_internal:count_internal+lphi(iii,jjj)-1) = integral_z(jjj); %phi_i  P controller
                elseif FLAG.PID(iii,jjj) == 2
                    phi_temp(:,count_internal:count_internal+lphi(iii,jjj)-1) = [z_buffer(jjj,1) integral_z(jjj) ]; %phi_i  PI controller
                elseif FLAG.PID(iii,jjj) == 3
                    phi_temp(:,count_internal:count_internal+lphi(iii,jjj)-1) = [z_buffer(jjj,1) integral_z(jjj) (z_buffer(jjj,1)-z_buffer(jjj,2))]; %phi_i  PID controller
                    %phi_temp(:,count_internal:count_internal+lphi(iii,jjj)-1) =  integral_z(jjj);
                end
            count_internal = count_internal + lphi(iii,jjj);
        end
        PHI(iii,count:count+size(phi_temp,2)-1) = phi_temp;
        count = count + size(phi_temp,2);
    end
    %Construct PHI

    u_bar(lu+1:end,:)       = u_bar(1:end-lu,:);
    u_bar(1:lu,:)           = u_in;

    PHI_bar(lu+1:end,:)     = PHI_bar(1:end-lu,:);
    PHI_bar(1:lu,:)         = PHI;
    
    %PHI_filt                = [zeros(lz,2*lu) filt_Nf]*PHI_bar;
    %u_filt                  = [zeros(lz,lu) filt_Nf]*u_bar;
    PHI_filt                = [zeros(lz,lu) filt_Nf]*PHI_bar;
    u_filt                  = filt_Nf*u_bar;
    
    X = [PHI_filt;PHI];
    gamma_inv = Rbar-Rbar*(X/(inv(P) + X'*Rbar*X))*X'*Rbar;
    Y = [ (PHI_filt*theta_buffer(:,1) + z_in - u_filt);(PHI*theta_buffer(:,1))];
    
    if FREEZE
         theta_out = theta_buffer(:,1);
    else  
         P = lambda_inv*P - lambda_inv*P*X'*gamma_inv*X*P;
         theta_out = theta_buffer(:,1) - P*X'*Rbar*Y ;
    end
    
    theta_buffer(:,2) = theta_buffer(:,1);
    theta_buffer(:,1) = theta_out;
    u_out             = PHI*theta_out; 
    P_out             = P;
    theta_id_out      = ThetaID;
    P_id_out          = P_id;
    lambda_id_out     = 1/betaID;

end
