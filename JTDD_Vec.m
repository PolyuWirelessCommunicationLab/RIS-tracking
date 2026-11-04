function [out] = JTDD_Vec(Y, param)
% Joint multiuser tracking and signal detection algorithm
% estimate AOA, delay, and Doppler
% INPUTS:
%   Y  - received signal at the BS size (M_BS, V_I, N_Fre, N_OFDM_symbol, BS_num)
%   K  - number of users

    %% parameter load
    flag_ini    = 1;
    while(flag_ini)
        ca_dim  = size(Y);  
        if length(ca_dim) == 4
            ca_dim(5) = 1;
        end
        dim_y   = prod(ca_dim(1:4));
        y       = reshape(Y,dim_y,ca_dim(5));
        T_out   = 20;
        T_VBI   = 300;
        
        bv_steer    = @(M_idx,pv) exp(1j*(M_idx-1).'*pv);
        f_VM_mk     = @(theta,mu,kappa) 1/(2*pi*besseli(0,kappa,1))*exp(kappa*cos(theta-mu));
        f_logVM_mk  = @(theta,mu,kappa) kappa*cos(theta-mu) - log(2*pi*besseli(0,kappa,1));
    
        c0          = 3e8;
        K_UE        = param.K;
        p_BS        = param.p_BS;
        p_RIS       = param.p_RIS;
        Pn          = param.Pn;
        M_RIS       = param.M_RIS;
        ss          = param.ss;
        sw          = param.sw;
        dt          = param.dt;
        th_IB       = param.theta_IB;
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
        A_true      = param.A_true;
        x           = param.x;

        ydd         = param.ydd;
        xdd0        = param.xdd;

        K_link      = K_UE*(num_BS + num_BS*num_RIS);
        K_link_UB   = K_UE*num_BS;
        K_link_UIB  = K_UE*num_BS*num_RIS;
    
        % J   = zeros([K_link,K_link]);
        % h   = zeros([K_link,1]);
        % w   = zeros([K_link,1]);
        % C   = zeros([K_link,K_link]);
        J   = zeros([K_UE*(1+num_RIS),K_UE*(1+num_RIS),num_BS]);
        h   = zeros([K_UE*(1+num_RIS),num_BS]);
        w   = zeros([K_UE*(1+num_RIS),num_BS]);
        C   = zeros([K_UE*(1+num_RIS),K_UE*(1+num_RIS),num_BS]);
        Kt  = zeros([T_VBI,num_BS]);
        mseB    = zeros([T_VBI,num_BS]);
        mse     = zeros([T_VBI,1]);
        
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
        pci_ftv_rt  = zeros([2,2,length(p_BS)]); % position prior inverse covariance
        vci_ftv_rt  = zeros([2,2,length(p_BS)]); % velocity prior inverse covariance
        for k = 1:K_UE
            pci_ftv_rt(:,:,k)   = eye(2)/pco_ftv_rt(:,:,k,1,1);
            vci_ftv_rt(:,:,k)   = eye(2)/vco_ftv_rt(:,:,k,1,1);
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
        mu_rt       = zeros([K_UE,4,1+num_RIS,num_BS,T_out]);  % K_UE x dim_resource (AOA(BS), AOA(IRS), Delay, Doppler) x links x Iteration
        kp_rt       = zeros([K_UE,4,1+num_RIS,num_BS,T_out]);
        eta_rt      = zeros([K_UE,4,1+num_RIS,num_BS,T_out]);
        ts_rt       = zeros([K_UE*(1+num_RIS),num_BS]);
        tsCo_rt     = P_link.^2;
        loc_error   = zeros([T_out+1,1]);
        vel_error   = zeros([T_out+1,1]);
        loc_error(1)= norm(p_Fu_bp-p_true,'fro'); 
        vel_error(1)= norm(vel_pri-v_true,'fro');
        mse_w       = zeros([T_out+1,1]);
        mse_w(1)    = norm(w_true,'fro');
        nmse_w      = zeros([T_out+1,1]);
        nmse_w(1)   = 0;
        mse_s       = zeros([T_out+1,1]);
        mse_s(1)    = norm(x,'fro');
        % st_pri      = true([K_UE*(1+num_RIS),num_BS]);
        st_pri      = param.st_pri;
        w_rec       = zeros([K_UE*(1+num_RIS),num_BS,T_out+1]);     w_rec(:,:,1) = w_true;
        mse_rec     = zeros([T_out+1,1]);

        st_pri_rec  = st_pri;

        pCi_Fu_to   = zeros([2,2,K_UE]);
        pCo_Fu_to   = zeros([2,2,K_UE]);
        p_Fu_to     = zeros([K_UE,2]);
        v_Fu_to     = zeros([K_UE,2]);
        vCo_Fu_to   = zeros([2,2,K_UE]);
        vCi_Fu_to   = zeros([2,2,K_UE]);

        flag_ini_param = 0;
    end
    

    % p_rec           = zeros([K_UE,2,num_BS+num_RIS,t_out]);
    t00 = 1;
    while(cond_out)

        %% Steering matrix update
        flag_stm = 1;
        while(flag_stm)
            for k = 1:K_UE
                for ib = 1:num_BS
                    param_rt(k,1,1,ib)      = norm(pm_ftv_rt(k,:,1,ib)-p_BS(ib,:), 2);   % distance for BS link
                    u_rt(k,:,1,ib)          = (pm_ftv_rt(k,:,1,ib)-p_BS(ib,:))/param_rt(k,1,1,ib); % direction vector for BS ib
                    tmpv_UB                 = [u_rt(k,2,1,ib), -u_rt(k,1,1,ib)];        % tangent vector for BS ib
                    v_rt(k,:,1,ib)          = tmpv_UB / norm(tmpv_UB,2);
                    param_rt(k,2,1,ib)      = u_rt(k,:,1,ib) * (ev_BS.');      
                    mu_rt(k,1,1,ib,t_out)   = pi * param_rt(k,2,1,ib);                % VM param for AOA of BS link
                    mod_num                 = floor((mu_rt(k,1,1,ib,t_out)+pi)/2/pi); % feasible region
                    mu_rt(k,1,1,ib,t_out)   = mu_rt(k,1,1,ib,t_out) - 2*pi*mod_num;
                    kp_rt(k,1,1,ib,t_out)   = (param_rt(k,1,1,ib)^2) / ...
                                            ((pi^2)*(1-param_rt(k,2,1,ib)^2)*(v_rt(k,:,1,ib)*pco_ftv_rt(:,:,k,1,ib)*(v_rt(k,:,1,ib).')));
                    eta_rt(k,1,1,ib,t_out)  = kp_rt(k,1,1,ib,t_out)*exp(1j*mu_rt(k,1,1,ib,t_out));
        
                    param_rt(k,3,1,ib)      = param_rt(k,1,1,ib)/c0;
                    mu_rt(k,3,1,ib,t_out)   = -2*pi*ss*param_rt(k,3,1,ib);
                    mu_rt(k,3,1,ib,t_out)   = mu_rt(k,3,1,ib,t_out); %- 2*pi*mod_num;
                    kp_rt(k,3,1,ib,t_out)   = c0^2/(4*(pi^2)*(ss^2)*(u_rt(k,:,1,ib)*pco_ftv_rt(:,:,k,1,ib)*(u_rt(k,:,1,ib).')));
                    eta_rt(k,3,1,ib,t_out)  = kp_rt(k,3,1,ib,t_out)*exp(1j*mu_rt(k,3,1,ib,t_out));
            
                    param_rt(k,4,1,ib)      = Nofs*vel_ftv_rt(k,:,1,ib)*(u_rt(k,:,1,ib).')/(lambda);
                    mu_rt(k,4,1,ib,t_out)   = -2*pi*dt*param_rt(k,4,1,ib);
                    mu_rt(k,4,1,ib,t_out)   = mu_rt(k,4,1,ib,t_out); %- 2*pi*mod_num;
                    kp_rt(k,4,1,ib,t_out)   = 1/((2*pi*dt*Nofs)^2*(u_rt(k,:,1,ib)*vco_ftv_rt(:,:,k,1,ib)*(u_rt(k,:,1,ib).')));
                    eta_rt(k,4,1,ib,t_out)  = kp_rt(k,4,1,ib,t_out)*exp(1j*mu_rt(k,4,1,ib,t_out));
                end
                
                for ib = 1:num_BS
                    for i = 1:num_RIS
                        ir = 1 + i;
                        mu_rt(k,1,ir,ib,t_out)  = th_IB(i+(ib-1)*num_RIS);
                        param_rt(k,1,ir,ib)     = norm(pm_ftv_rt(k,:,ir,ib)-p_RIS(i,:), 2);
                        u_rt(k,:,ir,ib)         = (pm_ftv_rt(k,:,ir,ib)-p_RIS(i,:))/param_rt(k,1,ir,ib);
                        tmpv_UI                 = [u_rt(k,2,ir,ib), -u_rt(k,1,ir,ib)];
                        v_rt(k,:,ir,ib)         = tmpv_UI / norm(tmpv_UI,2);
                        param_rt(k,2,ir,ib)     = u_rt(k,:,ir,ib) * (ev_RIS.');
                        mu_rt(k,2,ir,ib,t_out)  = pi * param_rt(k,2,ir,ib);
                        mu_rt(k,2,ir,ib,t_out)  = mu_rt(k,2,ir,ib,t_out); %- 2*pi*mod_num;
                        kp_rt(k,2,ir,ib,t_out)  = (param_rt(k,1,ir,ib)^2) / ...
                                                ((pi^2)*(1-param_rt(k,2,ir,ib)^2)*(v_rt(k,:,ir,ib)*pco_ftv_rt(:,:,k,ir,ib)*(v_rt(k,:,ir,ib).')));
                        eta_rt(k,2,ir,ib,t_out) = kp_rt(k,2,ir,ib,t_out)*exp(1j*mu_rt(k,2,ir,ib,t_out));
                
                        param_rt(k,3,ir,ib)     = (param_rt(k,1,ir,ib) + d_IB(i+(ib-1)*num_RIS)) / c0;
                        mu_rt(k,3,ir,ib,t_out)  = -2*pi*ss*param_rt(k,3,ir,ib);
                        kp_rt(k,3,ir,ib,t_out)  = c0^2/(4*(pi^2)*(ss^2)*(u_rt(k,:,ir,ib)*pco_ftv_rt(:,:,k,ir,ib)*(u_rt(k,:,ir,ib).')));
                        eta_rt(k,3,ir,ib,t_out) = kp_rt(k,3,ir,ib,t_out)*exp(1j*mu_rt(k,3,ir,ib,t_out));
                
                        param_rt(k,4,ir,ib)     = Nofs*vel_ftv_rt(k,:,ir,ib)*(u_rt(k,:,ir,ib).')/(lambda);
                        mu_rt(k,4,ir,ib,t_out)  = -2*pi*dt*param_rt(k,4,ir,ib);
                        mu_rt(k,4,ir,ib,t_out)     = mu_rt(k,4,ir,ib,t_out); %- 2*pi*mod_num;
                        kp_rt(k,4,ir,ib,t_out)  = 1/((2*pi*dt*Nofs)^2*(u_rt(k,:,ir,ib)*vco_ftv_rt(:,:,k,ir,ib)*(u_rt(k,:,ir,ib).')));
                        eta_rt(k,4,ir,ib,t_out) = kp_rt(k,4,ir,ib,t_out)*exp(1j*mu_rt(k,4,ir,ib,t_out));
                    end
                end
            end
            % eta_rt_sp   = zeros([K_UE,4,1+num_RIS,num_BS]);
            % for i = 1:num_BS
            %     eta_rt_sp(:,:,1,i)  = eta_rt(:,:,i,t_out);
            %     eta_rt_sp(:,:,2:(num_RIS+1),i)  = eta_rt(:,:,(num_BS+1):(num_BS+num_RIS),t_out);
            % end
        
            % Steering vector generation
            b_ud1       = zeros([ca_dim(1),K_link]);
            b_udS2      = zeros([M_RIS,K_link]);
            b_ud2       = zeros([ca_dim(2),K_link]);
            b_ud3       = zeros([ca_dim(3),K_link]);
            b_ud4       = zeros([ca_dim(4),K_link]);
            A           = zeros([prod(ca_dim),K_link]);
            mu_rt_BS    = zeros([4,K_UE*num_BS]);
            mu_rt_RIS   = zeros([4,K_UE*num_BS*num_RIS]);
            % Steering for BS link
            idx_BS              = 1:(K_UE*num_BS);
            num_BS_link         = K_UE*num_BS;
            mu_rt_BS(1,idx_BS)  = reshape(mu_rt(:,1,1,1:num_BS,t_out), [1,num_BS_link]);
            b_ud1(:,idx_BS)     = bv_steer(1:ca_dim(1), mu_rt_BS(1,idx_BS));
            b_ud2(:,idx_BS)     = ones([ca_dim(2),num_BS_link]);
            mu_rt_BS(3,idx_BS)  = reshape(mu_rt(:,3,1,1:num_BS,t_out), [1,num_BS_link]);
            b_ud3(:,idx_BS)     = bv_steer(1:ca_dim(3), mu_rt_BS(3,idx_BS));
            mu_rt_BS(4,idx_BS)  = reshape(mu_rt(:,4,1,1:num_BS,t_out), [1,num_BS_link]);
            b_ud4(:,idx_BS)     = bv_steer(1:ca_dim(4), mu_rt_BS(4,idx_BS));
            % Steering for RIS link
            idx_RIS             = (K_UE*num_BS+1):K_link;
            num_RIS_link        = K_UE*num_BS*num_RIS;
            idx_mu_RIS          = 1 + (1:num_RIS);
            b_ud1(:,idx_RIS)    = bv_steer(1:ca_dim(1), kron(th_IB,ones([1,K_UE])));
            mu_rt_RIS(2,:)      = reshape(mu_rt(:,2,idx_mu_RIS,:,t_out), [1,K_UE*num_RIS*num_BS]);
            b_udS2(:,idx_RIS)   = bv_steer(1:M_RIS, mu_rt_RIS(2,:));
            for ib = 1:num_BS
                for ir = 1:num_RIS
                    start               = (ib-1)*K_UE*num_RIS + (ir-1)*K_UE + 1;
                    ennd                = (ib-1)*K_UE*num_RIS + ir*K_UE;
                    idx_ibr             = K_UE*num_BS + (start:ennd);
                    b_ud2(:,idx_ibr)    = P(:,:,ir)*b_udS2(:,idx_ibr);
                end
            end
            % mu_rt_RIS(3,:)      = reshape(mu_rt(:,3,idx_mu_RIS,:,t_out), [1,K_UE*num_RIS*num_BS]) ... 
            %                     + kron(d_IB, ones([1,K_UE]))/c0*(-2*pi*ss);
            mu_rt_RIS(3,:)      = reshape(mu_rt(:,3,idx_mu_RIS,:,t_out), [1,K_UE*num_RIS*num_BS]);
            b_ud3(:,idx_RIS)    = bv_steer(1:ca_dim(3), mu_rt_RIS(3,:));
            mu_rt_RIS(4,:)      = reshape(mu_rt(:,4,idx_mu_RIS,:,t_out), [1,K_UE*num_RIS*num_BS]);
            b_ud4(:,idx_RIS)    = bv_steer(1:ca_dim(4), mu_rt_RIS(4,:));

            idx_pmt_BS  = zeros([1,K_link]);
            for ib = 1:num_BS
                id1             = ((ib-1)*K_UE*(1+num_RIS)+1):(ib*K_UE*(1+num_RIS));
                idx_pmt_BS(id1) = [((ib-1)*K_UE+1):(ib*K_UE), ((((ib-1)*K_UE*num_RIS+1):(ib*K_UE*num_RIS))+num_BS*K_UE)];
            end
            b_ud1   = b_ud1(:,idx_pmt_BS);
            b_udS2  = b_udS2(:,idx_pmt_BS);
            b_ud2   = b_ud2(:,idx_pmt_BS);
            b_ud3   = b_ud3(:,idx_pmt_BS);
            b_ud4   = b_ud4(:,idx_pmt_BS);

            b_ud1_kr    = kron( kron( kron(ones([ca_dim(4),1]),ones([ca_dim(3),1])), ones([ca_dim(2),1]) ), b_ud1 );
            b_ud2_kr    = kron( kron( kron(ones([ca_dim(4),1]),ones([ca_dim(3),1])), b_ud2 ), ones([ca_dim(1),1]) );
            b_ud3_kr    = kron( kron( kron(ones([ca_dim(4),1]),b_ud3), ones([ca_dim(2),1]) ), ones([ca_dim(1),1]) );
            b_ud4_kr    = kron( kron( kron(b_ud4,ones([ca_dim(3),1])), ones([ca_dim(2),1]) ), ones([ca_dim(1),1]) );
            A           = b_ud4_kr .* b_ud3_kr .* b_ud2_kr .* b_ud1_kr;
            A           = reshape(A, [dim_y,K_UE*(1+num_RIS),num_BS]);

            flag_stm = 0;
        end

        %% Start of Variational inference for parameters (q(theta) calculation)
        cont    = 1;
        t_VI    = 0;
        w_mean  = ts_rt;
        w_cov   = tsCo_rt;
    
        th_VBI = mu_rt(:,:,:,:,t_out);  
        kp_VBI = kp_rt(:,:,:,:,t_out);
        st_pri      = st_pri_rec;

        while(cont)
            t_VI    = t_VI + 1;
            s = st_pri;
            VirK    = sum(s,1);
            st_VBI  = s;
            K_UB    = sum(s(1:(K_UE),:),1);
            K_UI    = sum(s((K_UE+1):end,:),1);

            w       = zeros([K_UE*(1+num_RIS),num_BS]);
            wdd     = zeros([K_UE*(1+num_RIS),num_BS]);     
            wdd2    = zeros([K_UE*(1+num_RIS),num_BS]);
            xro     = zeros([dim_y,num_BS]);
            for ib = 1:num_BS
                A_use   = A(:,st_VBI(:,ib),ib);
                h(st_VBI(:,ib),ib)  = A_use'*(y(:,ib) - A_use*w_mean(st_VBI(:,ib),ib));
                tmpC                = (A_use'*A_use+diag(Pn./w_cov(st_VBI(:,ib),ib)));
                w(st_VBI(:,ib),ib)  = w_mean(st_VBI(:,ib),ib) + tmpC\h(st_VBI(:,ib),ib);
                C(st_VBI(:,ib),st_VBI(:,ib),ib)     = Pn*eye(VirK(ib))/tmpC;
                C(~st_VBI(:,ib),~st_VBI(:,ib),ib)   = diag(P_link(~st_VBI(:,ib),ib).^2);
                xro(:,ib)           = A(:,:,ib)*w(:,ib);

                wdd(st_VBI(:,ib),ib)  = tmpC\(A_use'*ydd(:,ib));
                A_use2  = A_true(:,st_VBI(:,ib),ib);
                wdd2(st_VBI(:,ib),ib)  = (A_use2'*A_use2+diag(Pn./w_cov(st_VBI(:,ib),ib)))\(A_use'*ydd(:,ib));
            end
            xdd     = wdd ./ P_link;
            xdd2    = wdd2 ./ P_link;

            C_kk = zeros([K_UE*(1+num_RIS),num_BS]);
            for ib = 1:num_BS
                C_kk(:,ib)  = real(diag(C(:,:,ib)));
                vomega      = 2*log(P_link(:,ib)) - log(C_kk(:,ib)) - (abs(w(:,ib)).^2)./C_kk(:,ib);
            
                s_small             = find(vomega > 3);
                st_pri_tmp          = true([K_UE*(1+num_RIS),1]);
                st_pri_tmp(s_small) = false([length(s_small),1]);
                st_pri(:,ib)        = st_pri(:,ib).*st_pri_tmp;
            end

    
            e_d_set = {};
            for d = 1:4
                e_d = zeros([ca_dim(d),1]); 
            end
            e_d = zeros([M_RIS,1]); e_d(2) = 1; e_d_set{5} = e_d;
            for i = 1:num_BS
                for k = 1:K_UE
                    if s(k,i) == true
                        inz     = (1:(K_UE*(1+num_RIS))).'; 
                        s_cj    = s(:,i);   s_cj(k) = false;
                        inzi    = inz(s_cj);
                        A_kb    = A(:,inzi,i);
                        tmp1    = 2/Pn * ( w(k,i)'*y(:,i) - A_kb*(C(inzi,k)+w(inzi,i)*(w(k,i)')) );
                        eta1    = tmp1;

                        th_rec  = th_VBI(k,:,1,i);  kp_rec = kp_VBI(k,:,1,i);
                        eta2_set{1} = eta_rt(k,1,1,i,t_out)*e_d_set{1};
                        eta2_set{2} = eta_rt(k,3,1,i,t_out)*e_d_set{3};
                        eta2_set{3} = eta_rt(k,4,1,i,t_out)*e_d_set{4};
                        [b_set, theta, kappa] = pntFreqVecEst(eta1, eta2_set, ca_dim, M_RIS, P, th_rec, kp_rec, 0);

                        th_VBI(k,:,1,i)     = theta;
                        kp_VBI(k,:,1,i)     = kappa;
                        l_idx               = (i-1)*(K_UE*(num_RIS+1)) + k;
                        b_ud1(:,l_idx)      = b_set{1};
                        b_ud2(:,l_idx)      = b_set{2};
                        b_ud3(:,l_idx)      = b_set{3};
                        b_ud4(:,l_idx)      = b_set{4};
                    end
                end
            end

            for ib = 1:num_BS
                for ir = 1:num_RIS
                    ki      = ir + 1;
                    pmt_dim = 1:4; pmt_dim(2)  = 1; pmt_dim(1)  = 2;
                    for k = 1:K_UE
                        idx_UEI         = k + ir*K_UE;
                        if s(idx_UEI,ib) ~= 0
                            inz     = (1:(K_UE*(1+num_RIS))).';
                            s_cj    = s(:,ib);   s_cj(idx_UEI) = false;
                            inzi    = inz(s_cj);
                            A_kb    = A(:,inzi,ib);
                            tmp1    = 2/Pn * ( w(idx_UEI,ib)'*y(:,ib) - A_kb*(C(inzi,idx_UEI,ib)+w(inzi,ib)*(w(idx_UEI,ib)')) );
                            eta1    = reshape(permute(reshape(tmp1,ca_dim(1:4)), pmt_dim), ca_dim(2), []);
                            eta2    = 1/Pn * (C(idx_UEI,idx_UEI,ib)+abs(w(idx_UEI,ib))^2) * (dim_y/ca_dim(2)) * vCP(:,ir);

                            th_rec  = th_VBI(k,:,ki,ib);  kp_rec = kp_VBI(k,:,ki,ib);
                            eta2_set{1} = eta_rt(k,2,ki,ib,t_out)*e_d_set{5} - eta2;
                            eta2_set{2} = eta_rt(k,3,ki,ib,t_out)*e_d_set{3};
                            eta2_set{3} = eta_rt(k,4,ki,ib,t_out)*e_d_set{4};
                            [b_set, theta, kappa] = pntFreqVecEst(eta1, eta2_set, ca_dim, M_RIS, P(:,:,ir), th_rec, kp_rec, 1);

                            th_VBI(k,:,ki,ib)   = theta;
                            kp_VBI(k,:,ki,ib)   = kappa;
                            idx_UI_BS           = k + ir*K_UE + (ib-1)*K_UE*(num_RIS+1);
                            b_ud1(:,idx_UI_BS)  = b_set{1};
                            b_ud2(:,idx_UI_BS)  = b_set{2};
                            b_ud3(:,idx_UI_BS)  = b_set{3};
                            b_ud4(:,idx_UI_BS)  = b_set{4};
                            b_udS2(:,idx_UI_BS) = b_set{5};
                        end
                    end
                end
            end
            
    

            b_ud1_kr    = kron( kron( kron(ones([ca_dim(4),1]),ones([ca_dim(3),1])), ones([ca_dim(2),1]) ), b_ud1 );
            b_ud2_kr    = kron( kron( kron(ones([ca_dim(4),1]),ones([ca_dim(3),1])), b_ud2 ), ones([ca_dim(1),1]) );
            b_ud3_kr    = kron( kron( kron(ones([ca_dim(4),1]),b_ud3), ones([ca_dim(2),1]) ), ones([ca_dim(1),1]) );
            b_ud4_kr    = kron( kron( kron(b_ud4,ones([ca_dim(3),1])), ones([ca_dim(2),1]) ), ones([ca_dim(1),1]) );
            A_tmp       = b_ud4_kr .* b_ud3_kr .* b_ud2_kr .* b_ud1_kr;
            A           = reshape(A_tmp, [dim_y,K_UE*(1+num_RIS),num_BS]);
    
    
            % Stopping criterion for the VBI process
            xr  = zeros([dim_y,num_BS]);
            for ib = 1:num_BS
                xr(:,ib)    = A(:,s(:,ib),ib)*w(s(:,ib),ib);
            end
            mseB(t_VI,:)  = 10*log10(sum(abs(xr-x_true).^2,1)./sum(abs(x_true).^2,1));
            mse(t_VI)   = 10*log10(norm(xr-x_true,'fro')^2/norm(x_true,'fro')^2);
            Kt(t_VI,:)    = VirK;
    
            if (norm(xr-xro)/norm(xro)<1e-8) || (norm(xro)==0&&norm(xr-xro)==0) || (t_VI >= T_VBI)
                cont = 0;
                mseB(t_VI:end,:)    = ones([size(mseB(t_VI:end,:),1),1])*mseB(t_VI,:);
                mse(t_VI:end)       = ones([size(mse(t_VI:end),1),1])   *mse(t_VI);
                % Kt(t_VI+1:end,:)    = Kt(t_VI,:);
            end
            xro = xr;
    
        end
        eta_VBI = kp_VBI .* exp(1j*th_VBI);

        mse_rec(t_out+1) = mse(end);
        if t00 > 0
            mse_iter_rec     = mse;
            t00 = 2;
        end

       
        %% extrinsic messages (theta) for fusion (ext(theta -> p_theta) = q(theta)/ext(p_theta -> theta))


        %% extrinsic messages (position) for fusion (ext(theta) => ext(position))


        %% position fusion (Gaussian) of BS link and RIS link (no prior information)


        %% extrinsic messages (velocity) for fusion


        %% transmit signal fusion


        t_out               = t_out + 1;
        loc_error(t_out)    = norm(p_ext_Fu_lt-p_true,'fro');
        vel_error(t_out)    = norm(vel_Fu-v_true,'fro');
        mse_w(t_out)        = norm(w-w_true,'fro');
        nmse_w(t_out)       = 10*log10(norm(w-w_true,'fro')^2/norm(w_true,'fro')^2);
        w_rec(:,:,t_out)    = w;
        mse_s(t_out)        = norm(ts_Fu-x,'fro').^2;
    
        xdd_re  = sum(reshape(xdd, [K_UE,(1+num_RIS)*num_BS]), 2);
        kt      = sum(reshape(st_VBI, [K_UE,(1+num_RIS)*num_BS]) , 2);
        xdd_est = xdd_re ./ kt;
        mse_xdd = norm(xdd_est-xdd0,'fro')^2;
        xdd2_re  = sum(reshape(xdd2, [K_UE,(1+num_RIS)*num_BS]), 2);
        kt      = sum(reshape(st_VBI, [K_UE,(1+num_RIS)*num_BS]) , 2);
        xdd2_est = xdd2_re ./ kt;
        mse_xdd2 = norm(xdd2_est-xdd0,'fro')^2;

        p_Fu_bp = p_ext_Fu_lt;
        v_Fu_bp = vel_Fu;
    
        if res < 1e-8 || t_out > T_out
            cond_out    = 0;
            loc_error(t_out:end)    = ones([size(loc_error(t_out:end),1),1])*loc_error(t_out);
            vel_error(t_out:end)    = ones([size(vel_error(t_out:end),1),1])*vel_error(t_out);
            mse_w(t_out:end)        = ones([size(mse_w(t_out:end),1),1])    *mse_w(t_out);
            nmse_w(t_out:end)       = ones([size(nmse_w(t_out:end),1),1])   *nmse_w(t_out);
            mse_s(t_out:end)        = ones([size(mse_s(t_out:end),1),1])    *mse_s(t_out);
            mse_rec(t_out:end)      = ones([size(mse_rec(t_out:end),1),1])  *mse_rec(t_out);
        end


        loc_error_pri(t_out)    = norm(p_ext_Fu_lt-pm_ftv_pri,'fro');
        vel_error_pri(t_out)    = norm(vel_Fu-vel_pri,'fro');
        
        p_Fu_to(:,:,t_out)      = p_ext_Fu_lt;
        pCo_Fu_to(:,:,:,t_out)  = pCo_ext_Fu_lt;
        pCi_Fu_to(:,:,:,t_out)  = pCi_ext_Fu_lt;
        v_Fu_to(:,:,t_out)      = vel_Fu;
        vCo_Fu_to(:,:,:,t_out)  = vCo_Fu;
        vCi_Fu_to(:,:,:,t_out)  = vCi_Fu;
        

    end
    

    out.pmean_pri       = p_ext_Fu_lt;
    out.pcov_pri        = pCo_ext_Fu_lt;
    out.pci_pri         = pCi_ext_Fu_lt;

    out.vel_pri         = vel_Fu;
    out.velcov_pri      = vCo_Fu;

    out.loc_error       = loc_error;
    out.vel_error       = vel_error;
    out.loc_error_fin   = loc_error(end);
    out.vel_error_fin   = vel_error(end);

    out.mse_s           = mse_s;
    out.mse_rec         = mse_rec;
    out.mse_iter_rec    = mse_iter_rec;
    out.mse_xdd         = mse_xdd;
    out.mse_xdd2        = mse_xdd2;

    p_st = 1;
    if any(find((st_pri-sw) < -0.5),'all')
        p_st = 0;
    end
    out.p_st = p_st;
    out.st   = st_pri;

end

function [f_value] = f_VMG(st_set, vel, vel_pri, vci_pri, th_set, kp_set, e_set, dt, Nofs, lambda)

    f_value = sum( (st_set.*(kp_set)) .* cos(-2*pi*dt*Nofs/lambda*vel*e_set-th_set), 2);
 
    f_value = f_value - 0.5*(vel-vel_pri)*vci_pri*(vel-vel_pri).';

end

function [vel, vCi, vCo] = pntVelFuEst(st_VBI_set, e_UBI_set, v_opt, vel_pri, vci_pri, th_set, kp_set, dt, Nofs, lambda)

    st_UE   = any(st_VBI_set, 'all');
    num_eta = size(th_set,2);
    sin_vel = zeros([1,num_eta]);
    cos_vel = zeros([1,num_eta]);
    
    %%%========Velocity is updated by eq. (26), (27)===============

    vel = v_opt;

end


function [b, theta, kappa] = pntFreqVecEst( eta1, eta2, ca_dim, M_RIS, P, th, kp, is_RIS )
% point estimation of the frequency vector (Vector form Gaussian Newton)

    theta   = th;
    kappa   = kp;
    
    %%%======Calculate the gradient and Hessian================
    %%%======Using eq. (46)==============


end


function [ k ] = Ainv( R )
% Returns the approximate solution of the equation R = A(k),
% where A(k) = I_1(k)/I_0(k) is the ration of modified Bessel functions of
% the first kind of first and zero order
% Uses the approximation from
%       Mardia & Jupp - Directional Statistics, Wiley 2000, pp. 85-86.
%
% When input R is a vector, the output is a vector containing the
% corresponding entries

    k   = R; % define A with same dimensions
    in1 = (R<.53); % indices of the entries < .53
    in3 = (R>=.85);% indices of the entries >= .85
    in2 = logical(1-in1-in3); % indices of the entries >=.53 and <.85
    R1  = R(in1); % entries < .53
    R2  = R(in2); % entries >=.53 and <.85
    R3  = R(in3); % entries >= .85
    
    % compute for the entries which are < .53
    if ~isempty(R1)
        t      = R1.*R1;
        k(in1) = R1 .* ( 2 + t + 5/6*t.*t );
    end
    % compute for the entries which are >=.53 and <.85
    if ~isempty(R2)
        k(in2) = -.4 + 1.39*R2 + 0.43./(1-R2);
    end
    % compute for the entries which are >= .85
    if ~isempty(R3)
        k(in3) = 1./( R3.*(R3-1).*(R3-3) );
    end

end
