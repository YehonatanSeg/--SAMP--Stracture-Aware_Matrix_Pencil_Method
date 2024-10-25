function [Poles, Left_modes,  MP_Amplitudes, NEW_Amplitudes,...
    SDD, GAP, EFF, KMEANS, MDL, AIC, CRB, params] = gen_data_SNR(noise_dist, Nreps, damps)
 
%% Settings
params.RANK_TYPE = 'effective';
params.NOSIE_TYPE = noise_dist; 
params.Nrep = Nreps;
params.SNR_RANGE = 0:2:30;
FREQ_DIST = 1.5;
INTERVAL_LENGTH = 10;
params.dt = 1/10;
params.interval_length = INTERVAL_LENGTH; 
N = 1 +params.interval_length / params.dt;
w1 = 2; 
params.true_freqs  = [ w1 + FREQ_DIST*pi/(params.dt*N), w1] ; 
params.true_amps   = [1, 1];
params.true_phases = [0, 0];

params.true_damps  = damps;
params.true_poles = exp(-params.true_damps*params.dt + 1i*params.true_freqs*params.dt);

Left_modes      = cell(length(params.SNR_RANGE), params.Nrep);
Poles           = cell(length(params.SNR_RANGE), params.Nrep);
MP_Amplitudes   = cell(length(params.SNR_RANGE), params.Nrep);
NEW_Amplitudes   = cell(length(params.SNR_RANGE), params.Nrep);

ind = 1;
for SNR = params.SNR_RANGE    
    
    % CRB
    SNR_linear = 10^(SNR/ 10);
    noise_var = sum(params.true_amps.^2) / SNR_linear;
    curr_CRB = CRB_unbiased_full_model( params.true_amps(1), params.true_amps(2), ...
                                        params.true_phases(1), params.true_phases(2), ...
                                        params.true_damps(1), params.true_damps(2), ...
                                        params.true_freqs(1), params.true_freqs(2), ...
                                        noise_var, N, params.dt ...
                                        );

    CRB.b1(ind)     = curr_CRB.b1;     
    CRB.b2(ind)     = curr_CRB.b2;    
    CRB.phi1(ind)   = curr_CRB.phi1;  
    CRB.phi2(ind)   = curr_CRB.phi2;
    CRB.alpha1(ind) = curr_CRB.alpha1; 
    CRB.alpha2(ind) = curr_CRB.alpha2;
    CRB.theta1(ind) = curr_CRB.theta1;
    CRB.theta2(ind) = curr_CRB.theta2;

    for rep = 1:params.Nrep
        disp(['Dampings = [', num2str(damps), ']', ' freq distance = ', num2str(FREQ_DIST), ' interval length = ', num2str(INTERVAL_LENGTH), ' SNR = ', num2str(SNR), ' rep = ', num2str(rep)])
        [params.signal, params.noise, params.noise_variance]  = generate_sum_of_k_complex_exp_with_add_noise(params.true_freqs, params.true_damps, params.true_phases, params.true_amps, params.dt, params.interval_length, SNR, params.NOSIE_TYPE);
      
        params.pencil_parameter = round(length(params.signal) /3);
        H = create_hankel_matrix(params.signal,params.pencil_parameter);
        [left_mode, mp_eigenvalues, b_MP, b_new] = get_MP_components(H, params.RANK_TYPE, []);

        % Full MP information for our methods
        Left_modes{ind,rep}      = left_mode;
        Poles{ind,rep}           = mp_eigenvalues;
        MP_Amplitudes{ind,rep}   = b_MP;
        NEW_Amplitudes{ind,rep}   = b_new;

        % Traditional MP methods ------------------------------------------
        [~,SDD_poles, b_SDD,   ~]        = get_MP_components(H, 'significant',[]);
        [~,gap_poles, b_gap, ~]        = get_MP_components(H, 'gap',[]);
        [~,eff_poles, b_eff, ~]        = get_MP_components(H, 'effective',[]);
        [~,kmeans_poles, b_kmeans, ~]  = get_MP_components(H, 'kmeans',[]);
        [~,mdl_poles, b_mdl, ~]        = get_MP_components(H, 'll_mdl', params);
        [~,aic_poles, b_aic, ~]        = get_MP_components(H, 'll_aic', params);

        SDD.Amplitudes.SV{ind,rep}     = b_SDD;
        SDD.Poles.SV{ind,rep}          = SDD_poles;
        SDD.NOF.SV(ind,rep)            = length(b_SDD);

        GAP.Amplitudes.SV{ind,rep}    = b_gap;
        GAP.Poles.SV{ind,rep}         = gap_poles;
        GAP.NOF.SV(ind,rep)           = length(b_gap);
        
        EFF.Amplitudes.SV{ind,rep}    = b_eff;
        EFF.Poles.SV{ind,rep}         = eff_poles;
        EFF.NOF.SV(ind,rep)           = length(b_eff);
            
        KMEANS.Amplitudes.SV{ind,rep} = b_kmeans;
        KMEANS.Poles.SV{ind,rep}      = kmeans_poles;
        KMEANS.NOF.SV(ind,rep)        = length(b_kmeans);

        MDL.Amplitudes.LL{ind,rep}    = b_mdl;
        MDL.Poles.LL{ind,rep}         = mdl_poles;
        MDL.NOF.LL(ind,rep)           = length(b_mdl);

        AIC.Amplitudes.LL{ind,rep}    = b_aic;
        AIC.Poles.LL{ind,rep}         = aic_poles;
        AIC.NOF.LL(ind,rep)           = length(b_aic);
        
    end
    
    ind = ind + 1;

end
