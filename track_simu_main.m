clc;
clear;

rng(133);


run_sys_param = 1;
while(run_sys_param)
    %% System parameters define
    M_BS        = 8;        % Antenna number 
    M_RIS       = 48;       % IRS element number
    dim_P       = 6;

    Q           = 12;       % Number of OFDM symbols
    L           = 4;        % Number of OFDM groups (each group repetive signal trans)
    LB          = 8;        % Number of repeative blocks for Doppler estimation
    QL          = Q*L;
    N_OFDM      = 3000;     % OFDM symbols in each coherence block
    TrackT      = 201;      % Number of coherence blocks'
    N_sdpl      = 10;       % Symbol number for Doppler estimation
    
    K_far       = 1;        % Number of far-field targets
    K_near      = 0;        % Number of near-field targets
    K           = K_far + K_near;
    
    N1          = 12;       % length of each subcarrier subblock
    N2          = 20;       % number of subcarrier subblocks
    N           = N1;       % Number of subcarrier
    N_CP        = N1;       % CP length for ISI mitigation
    BW_MHz      = 10;
    BW          = BW_MHz*(10^6); % bandwidth
    delta_f     = BW/N;
    dt          = 1/BW*(N_CP+N);
    dta_T       = 0.02;
    N_ofdm_syb  = 200;
    
    f_c         = 3e9;      % 3G common
    c0          = 3e8;
    lambda      = c0/f_c; 
    d_B         = lambda/2; % antenna space for BS antenna
    d_I         = lambda/2; % antenna space for IRS antenna

    Tx_power    = -50;
    Tx_apt      = 10^(Tx_power/10);
    sigma_z     = 10^(-170/10)*BW;
    Pn          = sigma_z/Tx_apt;
    
    num_BS  = 2;
    num_RIS = 2;
    p_BS    = zeros([num_BS,2]);
    p_RIS   = zeros([num_RIS,2]);
    for i = 1:num_BS
        if i == 1
            p_BS(i,:)   = [0,0];
        elseif i == 2
            p_BS(i,:)   = [90,0];
        end
    end
    for i = 1:num_RIS
        if i == 1
            p_RIS(i,:)  = [20,40];
        elseif i == 2
            p_RIS(i,:)  = [70,40];
        end
    end

    d_IB    = sqrt(vec(sum((repmat(p_RIS.',1,1,num_BS) ...
            - repmat(reshape(p_BS.',[2,1,num_BS]),1,num_RIS,1)).^2,1))).';
    th_IB   = zeros([1,num_RIS*num_BS]);
    e_BS    = [0;-1];
    e_RIS   = [1;0];
    for i = 1:num_BS
        for j = 1:num_RIS
            u_IB_tmp    = vec((p_RIS(j,:)-p_BS(i,:))/norm(p_RIS(j,:)-p_BS(i,:),2)).';
            th_IB(j+(i-1)*num_RIS)  = pi*u_IB_tmp*e_BS;
        end
    end
    P       = exp(1j*2*pi*rand([dim_P,M_RIS,num_RIS]));

    s_BI    = true([num_RIS,num_BS]);  
    s_BI(1,2) = false; s_BI(2,1) = false; %s_BI = false([num_RIS,num_BS]);
    sys_param.s_BI      = s_BI;

    %% System parameters import
    sys_param.M_BS      = M_BS;    % Antenna number 
    sys_param.M_RIS     = M_RIS;    % IRS element number
    
    sys_param.Q         = Q;    % Number of OFDM symbols
    sys_param.L         = L;   % Number of OFDM groups (each group repetive signal trans)
    sys_param.LB        = LB;    % Number of repeative blocks for Doppler estimation
    sys_param.QL        = QL;
    sys_param.N_OFDM    = N_OFDM; % OFDM symbols in each coherence block
    sys_param.TrackT    = TrackT;   % Number of coherence blocks
    sys_param.N_sdpl    = N_sdpl;
    
    sys_param.K_far     = K_far;    % Number of far-field targets
    sys_param.K_near    = K_near;    % Number of near-field targets
    sys_param.K         = K;
    
    sys_param.N1        = N1;   % length of each subcarrier subblock
    sys_param.N2        = N2;   % number of subcarrier subblocks
    sys_param.N         = N;    % Number of subcarrier
    sys_param.N_CP      = N_CP;  % CP length for ISI mitigation
    sys_param.BW_MHz    = BW_MHz;
    sys_param.BW        = BW; % bandwidth
    sys_param.delta_f   = delta_f;
    sys_param.dta_T     = dta_T;
    sys_param.dt        = dt;
    sys_param.N_ofdm_syb    = N_ofdm_syb;
    sys_param.Pn        = Pn;
    
    sys_param.f_c       = f_c; 
    sys_param.c0        = c0;
    sys_param.lambda    = lambda; 
    sys_param.d_B       = d_B; % antenna space for BS antenna
    sys_param.d_I       = d_I; % antenna space for IRS antenna

    sys_param.Tx_power      = Tx_power;
    sys_param.Tx_apt        = Tx_apt;
    sys_param.sigma_z       = sigma_z;
    param_text.Pn           = Pn;
    
    sys_param.num_BS        = num_BS;
    sys_param.num_RIS       = num_RIS;
    sys_param.p_BS          = p_BS;
    sys_param.p_RIS         = p_RIS;
    sys_param.d_IB          = d_IB;
    sys_param.theta_IB      = th_IB;
    sys_param.P             = P;    

    run_sys_param   = 0;
end


run_tga_lv = 1;
while(run_tga_lv)
    x_side = 50; y_side =  50;
    x_left = 20; y_down = -30;
    loc_target_start    = zeros([K_far,2]);
    vlc_target          = zeros([K_far,2]);
    acv_target          = zeros([K_far,2,TrackT]);
    vel_cov_ukn         = zeros([2,2,K_far]);
    
    %%%====generate velocity evolution in the trajectory=================

    %%%

    vlc_target_T    = zeros([K_far,2,TrackT]);
    loc_target_T    = zeros([K_far,2,TrackT]);
    for iT = 1:TrackT
        if iT == 1
            loc_target_T(:,:,iT)    = loc_target_start;
            vlc_target_T(:,:,iT)    = vlc_target;
        else
            for iK = 1:K_far
                vlc_target_T(iK,:,iT)   = vlc_target_T(iK,:,iT-1) + dta_T*acv_target(iK,:,iT-1);
                loc_target_T(iK,:,iT)   = loc_target_T(iK,:,iT-1) + dta_T*vlc_target_T(iK,:,iT-1);     
            end

        end

    end

    sys_param.loc_target_T  = loc_target_T;
    sys_param.vlc_target_T  = vlc_target_T;
    
    loc_target_T_pmt    = permute(loc_target_T, [2,3,1]);
    vlc_target_T_pmt    = permute(vlc_target_T, [2,3,1]);

    min_tau_UB  = 20/c0;
    max_tau_UB  = 77/c0;
    min_tau_UIB = 70/c0;
    max_tau_UIB = 131/c0;
    max_dpl     = 1/lambda*40;

    sys_param.min_tau_UB    = min_tau_UB;
    sys_param.max_tau_UB    = max_tau_UB;
    sys_param.min_tau_UIB   = min_tau_UIB;
    sys_param.max_tau_UIB   = max_tau_UIB;
    sys_param.max_dpl       = max_dpl;

    run_tga_lv   = 0;
end

bv_steer    = @(M_idx,pv) exp(1j*(M_idx-1).'*pv);


for iT = 1:201

    p_true  = loc_target_T(:,:,iT);
    param0  = zeros([K,4,2]);
    PL      = zeros([K*(1+num_RIS),num_BS]);
    u_BS    = zeros([K,2,num_BS]);
    u_RIS   = zeros([K,2,num_RIS]);
    Nss     = N1;       sys_param.Nss = Nss;
    ss      = BW/Nss;   sys_param.ss  = ss;
    dim_all = M_BS*dim_P*Nss*N_sdpl;

    for k = 1:K
        for i = 1:num_BS
            param0(k,1,i)   = norm(p_true(k,:)-p_BS(i,:),2);
            u_BS(k,:,i)     = (p_true(k,:)-p_BS(i,:))/param0(k,1,i);
            param0(k,2,i)   = pi*u_BS(k,:,i)*e_BS;
            param0(k,3,i)   = -2*pi*ss*param0(k,1,i)/c0;
            param0(k,4,i)   = -2*pi*dt*N_ofdm_syb/lambda*u_BS(k,:,i)*(vlc_target_T(k,:,iT)).';
        end
        for ki = 1:num_RIS
            i  = ki + num_BS;
            param0(k,1,i) = norm(p_true(k,:)-p_RIS(ki,:),2);
            u_RIS(k,:,ki) = (p_true(k,:)-p_RIS(ki,:))/param0(k,1,i);
            param0(k,2,i) = pi*u_RIS(k,:,ki)*e_RIS;
            param0(k,3,i) = -2*pi*ss*(param0(k,1,i))/c0;
            param0(k,4,i) = -2*pi*dt*N_ofdm_syb/lambda*u_RIS(k,:,ki)*(vlc_target_T(k,:,iT)).';
        end    
        for i = 1:num_BS
            PL(k,i) = 1/param0(k,1,i);
    
            for j = 1:num_RIS
                PL(k+j*K,i)  = 1/(param0(k,1,j+num_BS)*d_IB(j+(i-1)*num_RIS));
            end
        end
    end

    sys_param.param0   = param0;
    sys_param.v_true   = vlc_target_T(:,:,iT);
    sys_param.p_true   = loc_target_T(:,:,iT);

    x           = (randn([K,1]) + 1j*randn([K,1]))/sqrt(2);
    sys_param.x = x;   
    w       = kron(ones([1+num_RIS,num_BS]),x).*PL;
    sw      = ones([K*(1+num_RIS),num_BS]); 
    st_pri  = true([K*(1+num_RIS),num_BS]);
    for i = 1:num_BS
        idx1    = 1:K;
        if i == 1
            sw(idx1,i)  = zeros([K,1]);
        elseif i == 2
            sw(idx1,i)  = ones([K,1]);
        end
        sw(idx1,i)      = zeros([K,1]);
        st_pri(idx1,i)  = false([K,1]);
        for j = 1:num_RIS
            if s_BI(j,i) == 0
                idx         = (j*K+1):((j+1)*K);
                sw(idx,i)   = zeros([K,1]);
                st_pri(idx,i) = false([K,1]);
            end
        end
    end

    sys_param.st_pri    = st_pri;
    w                   = sw.*w;
    sys_param.sw        = sw;
    sys_param.P_link    = PL;
    sys_param.w_true    = w;

    %%%%%%==================================================================
    if iT == 1
        sys_param.pmean_pri    = loc_target_T(:,:,iT);
        sys_param.pcov_pri     = repmat(0.01*eye(2), 1,1,K);
        sys_param.vel_pri      = vlc_target_T(:,:,iT);
        sys_param.velcov_pri   = repmat(0.1*eye(2), 1,1,K) + vel_cov_ukn;
        BIM_old                = zeros([4*K,4*K]);
        for k = 1:K
            idx = ((k-1)*4+1):((k-1)*4+4);
            BIM_old(idx,idx)    = inv([[0.01*eye(2),zeros([2,2])]; ...
                                   [zeros([2,2]),0.1*eye(2)+vel_cov_ukn(:,:,k)]]);
        end
    else
        if any(isnan(outJTDD.pmean_pri),'all') || any(isnan(outJTDD.loc_error),'all')
            sys_param.pmean_pri    = loc_target_T(:,:,iT-1);
            sys_param.pcov_pri     = repmat(1*eye(2), 1,1,K);
            sys_param.vel_pri      = vlc_target_T(:,:,iT-1);
            sys_param.velcov_pri   = repmat(5*eye(2), 1,1,K) + vel_cov_ukn;
        else
            sys_param.pmean_pri    = outJTDD.pmean_pri;
            sys_param.pcov_pri     = outJTDD.pcov_pri;
            sys_param.vel_pri      = outJTDD.vel_pri;
            sys_param.velcov_pri   = outJTDD.velcov_pri  + vel_cov_ukn;
        end
        BIM_old = BIM_sre(:,:,iT-1,iR);
    end
    %%%%%%==================================================================
    

    
    A       = zeros([dim_all,K*(1+num_RIS)*num_BS]);
    as1     = zeros([M_BS,   K*(1+num_RIS)*num_BS]);
    as2     = zeros([dim_P,  K*(1+num_RIS)*num_BS]);
    as3     = zeros([Nss,    K*(1+num_RIS)*num_BS]);
    as4     = zeros([N_sdpl, K*(1+num_RIS)*num_BS]);    
    for k = 1:K
        for i = 1:num_BS
            k_UB        = (i-1)*K*(1+num_RIS) + k;
            as1(:,k_UB) = bv_steer(1:M_BS,squeeze(param0(k,2,i)));
            as2(:,k_UB) = ones([dim_P,1]);
            as3(:,k_UB) = bv_steer(1:Nss,squeeze(param0(k,3,i)));
            as4(:,k_UB) = bv_steer(1:N_sdpl,squeeze(param0(k,4,i)));
            A(:,k_UB)   = kron(kron(kron(as4(:,k_UB),as3(:,k_UB)),as2(:,k_UB)),as1(:,k_UB));
    
            for j = 1:num_RIS
                k_UI        = (i-1)*K*(1+num_RIS) + j*K + k;
                as1(:,k_UI) = bv_steer(1:M_BS,th_IB(j+(i-1)*num_RIS));
                as2(:,k_UI) = P(:,:,j)*bv_steer(1:M_RIS,param0(k,2,j+num_BS));
                as3(:,k_UI) = bv_steer(1:Nss,param0(k,3,j+num_BS)-2*pi*ss/c0*d_IB(j+(i-1)*num_RIS));
                as4(:,k_UI) = bv_steer(1:N_sdpl,param0(k,4,j+num_BS));
                A(:,k_UI)   = kron(kron(kron(as4(:,k_UI),as3(:,k_UI)),as2(:,k_UI)),as1(:,k_UI));
            end
        end
    end
    A   = reshape(A, [dim_all,K*(1+num_RIS),num_BS]);
    
    sys_param.s_link0   = ones([K*(1+num_RIS),1]);
    y_nn                = zeros([dim_all,num_BS]);
    y                   = zeros([dim_all,num_BS]);
    noise               = sqrt(param_text.Pn)/sqrt(2)*(randn(dim_all,num_BS) + 1j*randn(dim_all,num_BS));
    for i = 1:num_BS
        y_nn(:,i)   = A(:,:,i)*w(:,i);
        y(:,i)      = y_nn(:,i) + noise(:,i);
    end
    Y                   = reshape(y, [M_BS,dim_P,Nss,N_sdpl,num_BS]);
    sys_param.x_true    = reshape(y_nn, [dim_all,num_BS]);
    sys_param.A_true    = A;


    xdd         = (randn([K,1]) + 1j*randn([K,1]))/sqrt(2);
    wdd         = kron(ones([1+num_RIS,num_BS]),xdd).*PL;
    y_nn_dd     = zeros([dim_all,num_BS]);
    y_dd        = zeros([dim_all,num_BS]);
    noise_dd    = sqrt(param_text.Pn)/sqrt(2)*(randn(dim_all,num_BS) + 1j*randn(dim_all,num_BS));
    for i = 1:num_BS
        y_nn_dd(:,i)   = A(:,:,i)*wdd(:,i);
        y_dd(:,i)      = y_nn_dd(:,i) + noise_dd(:,i);
    end
    sys_param.ydd   = y_dd;    
    sys_param.xdd   = xdd;

    [BCRB, BIM, BIM_state, BCRB_p,BCRB_v,BCRB_w] = BCRB(Y, sys_param, BIM_old);
    outJTDD1             = JTDD_Vec(Y, sys_param);
    p_st_T(:,iT)         = outJTDD.p_st;
    st_est_T(:,:,iT)     = outJTDD.st;
    loc_error_T(:,iT)    = outJTDD.loc_error;
    vel_error_T(:,iT)    = outJTDD.vel_error;
    loc_est_T(:,:,iT)    = outJTDD.pmean_pri;
    vel_est_T(:,:,iT)    = outJTDD.vel_pri;
    mseS_est_T(:,:,iT)   = outJTDD.mse_s;
    loc_BCRB_T(:,iT)     = BCRB_p;
    vel_BCRB_T(:,iT)     = BCRB_v;
    mseS_BCRB_T(:,iT)    = BCRB_w;
    BIM_sre(:,:,iT)      = BIM_state;
    mse_rec_T(:,iT)      = outJTDD.mse_rec;
    mse_rec_iter_T(:,iT) = outJTDD.mse_iter_rec;


end






