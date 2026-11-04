function [BCRB,BIM,BIM_state,BCRB_p,BCRB_v,BCRB_w] = BCRB(Y, param)
% BCRB: the BCRB calculation for user state and transmit signal

    %% parameter load
    flag_ini    = 1;
    while(flag_ini)
        ca_dim  = size(Y);  
        if length(ca_dim) == 4
            ca_dim(5) = 1;
        end
        dim_y   = prod(ca_dim(1:4));
        y       = reshape(Y,dim_y,ca_dim(5));
        
        bv_steer    = @(M_idx,pv) exp(1j*(M_idx-1).'*pv);
        f_VM_mk     = @(theta,mu,kappa) 1/(2*pi*besseli(0,kappa,1))*exp(kappa*cos(theta-mu));
        f_logVM_mk  = @(theta,mu,kappa) kappa*cos(theta-mu) - log(2*pi*besseli(0,kappa,1));
    
        c0          = 3e8;
        K_UE        = param.K;
        p_BS        = param.p_BS;
        p_RIS       = param.p_RIS;
        Pn          = param.Pn;
        th_IB       = param.theta_IB;
        M_RIS       = param.M_RIS;
        ss          = param.ss;
        sw          = param.sw;
        dt          = param.dt;
        d_IB        = param.d_IB;
        pmean_pri   = param.pmean_pri;
        pcov_pri    = param.pcov_pri;
        vel_pri     = param.vel_pri;
        velcov_pri  = param.velcov_pri;
        x_true      = param.x_true;
        dta_T       = param.dta_T;
        lambda      = param.lambda;
        P_link      = param.P_link;
        w_true      = param.w_true;
        Nofs        = param.N_ofdm_syb;
        p_true      = param.p_true;
        v_true      = param.v_true;
        num_RIS     = param.num_RIS;
        num_BS      = param.num_BS;
        s_link0     = param.s_link0;
        x           = param.x;

        K_link      = K_UE*(num_BS + num_BS*num_RIS);
        K_link_UB   = K_UE*num_BS;
        K_link_UIB  = K_UE*num_BS*num_RIS;
        
        P       = param.P;
        vCP     = zeros([M_RIS,num_RIS]);     % for theta(UE-RIS) update
        for iR = 1:num_RIS
            CP          = P(:,:,iR)'*P(:,:,iR);
            vCP(1,iR)   = ca_dim(2)*M_RIS;
            for i = 2:M_RIS
                for j = 1:(M_RIS-i+1)
                    vCP(i,iR) = vCP(i,iR) + 2*CP(j+i-1,j);
                end
            end
        end

        flag_ini = 0;
    end

    % Initialization: theta, tau, doppler from position and velocity
    flag_ini_param = 1;
    while(flag_ini_param)
        
        pm_ftv_rt   = pmean_pri + vel_pri*dta_T;    pm_ftv_pri  = pm_ftv_rt; %position prior mean     
        p_Fu_bp     = pm_ftv_rt;
        num_ftv     = num_BS + num_RIS;
        pm_ftv_rt   = repmat(pm_ftv_rt,  1,1,1+num_RIS,num_BS);
        pco_ftv_rt  = pcov_pri  + velcov_pri*dta_T; pco_ftv_pri = pco_ftv_rt;%position prior covariance
        pco_ftv_rt  = repmat(pco_ftv_rt, 1,1,1,1+num_RIS,num_BS);
        vel_ftv_rt  = vel_pri;    % velocity prior mean                              
        v_Fu_bp     = vel_ftv_rt;
        vel_ftv_rt  = repmat(vel_ftv_rt, 1,1,1+num_RIS,num_BS);
        velcov_rt   = velcov_pri;                   vco_ftv_rt  = repmat(velcov_rt,  1,1,1,1+num_RIS,num_BS); % velocity prior covariance
        pci_ftv_rt  = zeros([length(p_BS),length(p_BS)]); % position prior inverse covariance
        vci_ftv_rt  = zeros([length(p_BS),length(p_BS)]); % velocity prior inverse covariance
        for k = 1:K_UE
            pci_ftv_rt(:,:,k)   = eye(2)/pco_ftv_rt(:,:,k,1);
            vci_ftv_rt(:,:,k)   = eye(2)/vco_ftv_rt(:,:,k,1);
        end
        pci_pri     = pci_ftv_rt;
        vci_pri     = vci_ftv_rt;
        pci_ftv_rt  = repmat(pci_ftv_rt,  1,1,1,1+num_RIS,num_BS);
        vci_ftv_rt  = repmat(vci_ftv_rt,  1,1,1,1+num_RIS,num_BS);
        ev_BS = [0,-1];   ev_RIS = [1,0];
    
    
        cond_out    = 1;
        t_out       = 1;
        param_rt    = zeros([K_UE,4,1+num_RIS,num_BS]);  % distance x AOA(BS/IRS) x delay x doppler x links
        v_rt        = zeros([K_UE,2,1+num_RIS,num_BS]);  % K_UE x dim_cordinate x links
        u_rt        = zeros([K_UE,2,1+num_RIS,num_BS]);
        mu_rt       = zeros([K_UE,5,1+num_RIS,num_BS]);  % K_UE x dim_resource (AOA(BS), AOA(IRS), Delay, Doppler) x links x Iteration
        kp_rt       = zeros([K_UE,5,1+num_RIS,num_BS]);
        eta_rt      = zeros([K_UE,5,1+num_RIS,num_BS]);
        ts_rt       = zeros([K_UE*(1+num_RIS),num_BS]);
        tsCo_rt     = P_link.^2;

        pCi_Fu_to   = zeros([2,2,K_UE]);
        pCo_Fu_to   = zeros([2,2,K_UE]);
        p_Fu_to     = zeros([K_UE,2]);
        v_Fu_to     = zeros([K_UE,2]);
        vCo_Fu_to   = zeros([2,2,K_UE]);
        vCi_Fu_to   = zeros([2,2,K_UE]);

        flag_ini_param = 0;
    end

    bv_steer    = @(M_idx,pv) exp(1j*(M_idx-1).'*pv);
    bv_steer_th = @(M_idx,pv_sin) 1j*(M_idx-1).'*pv_sin;
    

    %% Steering matrix update
    flag_stm = 1;
    while(flag_stm)
        for k = 1:K_UE
            for ib = 1:num_BS
                param_rt(k,1,1,ib)      = norm(pm_ftv_rt(k,:,1,ib)-p_BS(ib,:), 2);   % distance for BS link
                u_rt(k,:,1,ib)          = (pm_ftv_rt(k,:,1,ib)-p_BS(ib,:))/param_rt(k,1,1,ib); % direction vector for BS ib
                % param_rt(k,1,1,ib)      = norm(p_true(k,:)-p_BS(ib,:), 2);   % distance for BS link
                % u_rt(k,:,1,ib)          = (p_true(k,:)-p_BS(ib,:))/param_rt(k,1,1,ib); % direction vector for BS ib
                tmpv_UB                 = [u_rt(k,2,1,ib), -u_rt(k,1,1,ib)];        % tangent vector for BS ib
                v_rt(k,:,1,ib)          = tmpv_UB / norm(tmpv_UB,2);
                param_rt(k,2,1,ib)      = u_rt(k,:,1,ib) * (ev_BS.');      
                mu_rt(k,1,1,ib)         = pi * param_rt(k,2,1,ib);                % VM param for AOA of BS link
                mu_rt(k,5,1,ib)         = pi * sqrt(1-param_rt(k,2,1,ib)^2);
                kp_rt(k,1,1,ib)         = (param_rt(k,1,1,ib)^2) / ...
                                            ((pi^2)*(1-param_rt(k,2,1,ib)^2)*(v_rt(k,:,1,ib)*pco_ftv_rt(:,:,k,1,ib)*(v_rt(k,:,1,ib).')));
                eta_rt(k,1,1,ib)        = kp_rt(k,1,1,ib)*exp(1j*mu_rt(k,1,1,ib));
        
                param_rt(k,3,1,ib)      = param_rt(k,1,1,ib)/c0;
                mu_rt(k,3,1,ib)         = -2*pi*ss*param_rt(k,3,1,ib);
                mu_rt(k,3,1,ib)         = mu_rt(k,3,1,ib); %- 2*pi*mod_num;
                kp_rt(k,3,1,ib)         = c0^2/(4*(pi^2)*(ss^2)*(u_rt(k,:,1,ib)*pco_ftv_rt(:,:,k,1,ib)*(u_rt(k,:,1,ib).')));
                eta_rt(k,3,1,ib)        = kp_rt(k,3,1,ib)*exp(1j*mu_rt(k,3,1,ib));
            
                param_rt(k,4,1,ib)      = vel_ftv_rt(k,:,1,ib)*(u_rt(k,:,1,ib).')/(lambda);
                % param_rt(k,4,1,ib)      = v_true(k,:)*(u_rt(k,:,1,ib).')/(lambda);
                mu_rt(k,4,1,ib)         = -2*pi*dt*Nofs*param_rt(k,4,1,ib);
                mu_rt(k,4,1,ib)         = mu_rt(k,4,1,ib); %- 2*pi*mod_num;
                kp_rt(k,4,1,ib)         = 1/((2*pi*dt*Nofs)^2*(u_rt(k,:,1,ib)*vco_ftv_rt(:,:,k,1,ib)*(u_rt(k,:,1,ib).')));
                eta_rt(k,4,1,ib)        = kp_rt(k,4,1,ib)*exp(1j*mu_rt(k,4,1,ib));
            end
                
            for ib = 1:num_BS
                for i = 1:num_RIS
                    ir = 1 + i;
                    param_rt(k,1,ir,ib)     = norm(pm_ftv_rt(k,:,ir,ib)-p_RIS(i,:), 2);
                    u_rt(k,:,ir,ib)         = (pm_ftv_rt(k,:,ir,ib)-p_RIS(i,:))/param_rt(k,1,ir,ib);
                    % param_rt(k,1,ir,ib)     = norm(p_true(k,:)-p_RIS(i,:), 2);
                    % u_rt(k,:,ir,ib)         = (p_true(k,:)-p_RIS(i,:))/param_rt(k,1,ir,ib);
                    tmpv_UI                 = [u_rt(k,2,ir,ib), -u_rt(k,1,ir,ib)];
                    v_rt(k,:,ir,ib)         = tmpv_UI / norm(tmpv_UI,2);
                    param_rt(k,2,ir,ib)     = u_rt(k,:,ir,ib) * (ev_RIS.');
                    mu_rt(k,2,ir,ib)        = pi * param_rt(k,2,ir,ib);
                    mu_rt(k,5,ir,ib)        = pi * sqrt(1-param_rt(k,2,ir,ib)^2);
                    kp_rt(k,2,ir,ib)        = (param_rt(k,1,ir,ib)^2) / ...
                                                ((pi^2)*(1-param_rt(k,2,ir,ib)^2)*(v_rt(k,:,ir,ib)*pco_ftv_rt(:,:,k,ir,ib)*(v_rt(k,:,ir,ib).')));
                    eta_rt(k,2,ir,ib)       = kp_rt(k,2,ir,ib)*exp(1j*mu_rt(k,2,ir,ib));
                
                    param_rt(k,3,ir,ib)     = (param_rt(k,1,ir,ib)) / c0;
                    mu_rt(k,3,ir,ib)        = -2*pi*ss*param_rt(k,3,ir,ib);
                    kp_rt(k,3,ir,ib)        = c0^2/(4*(pi^2)*(ss^2)*(u_rt(k,:,ir,ib)*pco_ftv_rt(:,:,k,ir,ib)*(u_rt(k,:,ir,ib).')));
                    eta_rt(k,3,ir,ib)       = kp_rt(k,3,ir,ib)*exp(1j*mu_rt(k,3,ir,ib));
                
                    param_rt(k,4,ir,ib)     = vel_ftv_rt(k,:,ir,ib)*(u_rt(k,:,ir,ib).')/(lambda);
                    % param_rt(k,4,ir,ib)     = v_true(k,:)*(u_rt(k,:,ir,ib).')/(lambda);
                    mu_rt(k,4,ir,ib)        = -2*pi*dt*Nofs*param_rt(k,4,ir,ib);
                    mu_rt(k,4,ir,ib)        = mu_rt(k,4,ir,ib,t_out); %- 2*pi*mod_num;
                    kp_rt(k,4,ir,ib)        = 1/((2*pi*dt*Nofs)^2*(u_rt(k,:,ir,ib)*vco_ftv_rt(:,:,k,ir,ib)*(u_rt(k,:,ir,ib).')));
                    eta_rt(k,4,ir,i)        = kp_rt(k,4,ir,ib)*exp(1j*mu_rt(k,4,ir,ib));
                end
            end
        end

        b_ud1_UB    = zeros([ca_dim(1),K_UE,num_BS]);
        b_ud2_UB    = zeros([ca_dim(2),K_UE,num_BS]);
        b_ud3_UB    = zeros([ca_dim(3),K_UE,num_BS]);
        b_ud4_UB    = zeros([ca_dim(4),K_UE,num_BS]);
        mu1_BS      = reshape(mu_rt(:,1,1,1:num_BS), [1,K_UE*num_BS]);
        b_ud1_UB    = reshape(bv_steer(1:ca_dim(1), mu1_BS), [],K_UE,num_BS);
        mu5_BS      = reshape(mu_rt(:,5,1,1:num_BS), [1,K_UE*num_BS]);
        db_ud1_UB1  = reshape( ((1:ca_dim(1)).'-1) * (-1j * mu5_BS), [],K_UE,num_BS);
        db_ud1_UB   = db_ud1_UB1 .* b_ud1_UB;
        b_ud2_UB    = ones([ca_dim(2),K_UE,num_BS]);
        mu3_BS      = reshape(mu_rt(:,3,1,1:num_BS), [1,K_UE*num_BS]);
        b_ud3_UB    = reshape(bv_steer(1:ca_dim(3), mu3_BS), [],K_UE,num_BS);
        db_ud3_UB1  = reshape( ((1:ca_dim(3)).'-1) * (1j * (-2*pi*ss)/c0 * ones([1,K_UE*num_BS])), [],K_UE,num_BS);
        db_ud3_UB   = db_ud3_UB1 .* b_ud3_UB;
        mu4_BS      = reshape(mu_rt(:,4,1,1:num_BS), [1,K_UE*num_BS]);
        b_ud4_UB    = reshape(bv_steer(1:ca_dim(4), mu4_BS), [],K_UE,num_BS);
        db_ud4_UB1  = reshape( ((1:ca_dim(4)).'-1) * (1j * (-2*pi*dt*Nofs) * ones([1,K_UE*num_BS])), [],K_UE,num_BS);
        db_ud4_UB   = db_ud4_UB1 .* b_ud4_UB;

        b_ud1_UIB   = zeros([ca_dim(1),K_UE,num_RIS,num_BS]);
        b_ud2_UIB   = zeros([ca_dim(2),K_UE,num_RIS,num_BS]);
        db_ud2_UIB  = zeros([ca_dim(2),K_UE,num_RIS,num_BS]);
        b_ud3_UIB   = zeros([ca_dim(3),K_UE,num_RIS,num_BS]);
        b_ud4_UIB   = zeros([ca_dim(4),K_UE,num_RIS,num_BS]);
        idx_ir      = 2:(num_RIS+1);
        mu1_RIS     = kron(th_IB,ones([1,K_UE]));
        b_ud1_UIB   = reshape(bv_steer(1:ca_dim(1), mu1_RIS), [],K_UE,num_RIS,num_BS);
        mu2_RIS     = reshape(mu_rt(:,2,idx_ir,1:num_BS), [1,K_UE*num_RIS*num_BS]);
        b_udS2_UIB  = reshape(bv_steer(1:M_RIS, mu2_RIS), [],K_UE,num_RIS,num_BS);
        mu5_RIS     = reshape(mu_rt(:,5,idx_ir,1:num_BS), [1,K_UE*num_RIS*num_BS]);
        db_udS2_UIB1= reshape( ((1:M_RIS).'-1) * (-1j * mu5_RIS), [],K_UE,num_RIS,num_BS);
        db_udS2_UIB = db_udS2_UIB1 .* b_udS2_UIB;
        for ib = 1:num_BS
            for ir = 1:num_RIS
                b_ud2_UIB(:,:,ir,ib)    = P(:,:,ir)*b_udS2_UIB(:,:,ir,ib);
                db_ud2_UIB(:,:,ir,ib)   = P(:,:,ir)*db_udS2_UIB(:,:,ir,ib);
            end
        end
        mu3_RIS     = reshape(mu_rt(:,3,idx_ir,1:num_BS), [1,K_UE*num_RIS*num_BS]) ...
                    + kron(d_IB, ones([1,K_UE]))/c0*(-2*pi*ss);
        b_ud3_UIB   = reshape(bv_steer(1:ca_dim(3), mu3_RIS), [],K_UE,num_RIS,num_BS);
        db_ud3_UIB1 = reshape( ((1:ca_dim(3)).'-1) * (1j * (-2*pi*ss)/c0 * ones([1,K_UE*num_RIS*num_BS])), [],K_UE,num_RIS,num_BS);
        db_ud3_UIB  = db_ud3_UIB1 .* b_ud3_UIB;
        mu4_RIS     = reshape(mu_rt(:,4,idx_ir,1:num_BS), [1,K_UE*num_RIS*num_BS]);
        b_ud4_UIB   = reshape(bv_steer(1:ca_dim(4), mu4_RIS), [],K_UE,num_RIS,num_BS);
        db_ud4_UIB1 = reshape( ((1:ca_dim(4)).'-1) * (1j * (-2*pi*dt*Nofs) * ones([1,K_UE*num_RIS*num_BS])), [],K_UE,num_RIS,num_BS);
        db_ud4_UIB  = db_ud4_UIB1 .* b_ud4_UIB;

        flag_stm = 0;
    end

    %%%======BIM calculation==========

    %%%======Jacobian calculation=======



    BIM         = B_mea_state + B_pri;
    BIM_state   = Jab*(2/Pn*real(B1))*(Jab.') + Bpv_pri_tp1;
    BIM_w       = 2/Pn*real(B3) + eye(K_UE);
    BCRB        = inv(BIM);
    BCRB_state  = inv(BIM_state);
    BCRB_w      = sum(diag(inv(BIM_w)));
    
    idx_p           = repmat([true;true;false;false], K_UE,1);
    idx_v           = repmat([false;false;true;true], K_UE,1);
    diag_BCRB_state = diag(BCRB_state);
    BCRB_p          = sqrt(sum(diag_BCRB_state(idx_p)));
    BCRB_v          = sqrt(sum(diag_BCRB_state(idx_v)));



end

