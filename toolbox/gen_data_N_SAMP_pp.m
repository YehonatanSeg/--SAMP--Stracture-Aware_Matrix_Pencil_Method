function [Poles, Left_modes, Amplitudes,...
    SDD, GAP, EFF, KMEANS, MDL, AIC, MAP, EVT, params] = gen_data_N_SAMP_pp(noise_distribution, Nreps, damps)
%% Settings
params.RANK_TYPE = 'effective';
params.NOISE_TYPE = noise_distribution;
params.Nrep = Nreps;
params.SNR_RANGE = 10;
params.FREQ_DIST = 1.5;
params.INTERVAL_LENGTH_RANGE = 2:2:20;
SNR = params.SNR_RANGE;
params.dt = 1/15;
params.true_amps   = [1, 1, 1, 1];
params.true_phases = [0, 0, 0, 0];
params.true_damps  = damps;
w1 = 2; 
w2 = 4;

Left_modes   = cell(length(params.INTERVAL_LENGTH_RANGE), params.Nrep);
Poles        = cell(length(params.INTERVAL_LENGTH_RANGE), params.Nrep);
Amplitudes   = cell(length(params.INTERVAL_LENGTH_RANGE), params.Nrep);

ind = 1;
for intLength = params.INTERVAL_LENGTH_RANGE 
    
    params.interval_length = intLength; 
    N = 1 +params.interval_length / params.dt;
    params.true_freqs  = [w2 + params.FREQ_DIST*pi/(params.dt*N), w2, w1 + params.FREQ_DIST*pi/(params.dt*N), w1]; 
    params.true_poles = exp(-params.true_damps*params.dt + 1i*params.true_freqs*params.dt);
    
    for rep = 1:params.Nrep
        disp(['Dampings = [', num2str(damps), ']', ' interval length = ', num2str(intLength), ' SNR = ', num2str(SNR), ' rep = ', num2str(rep)])

        [params.signal, params.noise, params.noise_variance]  = generate_sum_of_k_complex_exp_with_add_noise(params.true_freqs, params.true_damps, params.true_phases, params.true_amps, params.dt, params.interval_length, SNR, params.NOISE_TYPE);

        params.pencil_parameter = round(length(params.signal) /3);
        H = create_hankel_matrix(params.signal,params.pencil_parameter);
        [left_mode, mp_eigenvalues, ~, b_new] = get_MP_components(H, params.RANK_TYPE, []);
        
        % Full MP information for our methods
        Left_modes{ind,rep}   = left_mode;
        Poles{ind,rep}        = mp_eigenvalues;
        Amplitudes{ind,rep}   = b_new;

        % Traditional MP methods ------------------------------------------
        [~,SDD_poles, b_SDD,   ~]      = get_MP_components(H, 'significant',[]);
        [~,gap_poles, b_gap, ~]        = get_MP_components(H, 'gap',[]);
        [~,eff_poles, b_eff, ~]        = get_MP_components(H, 'effective',[]);
        [~,kmeans_poles, b_kmeans, ~]  = get_MP_components(H, 'kmeans',[]);
        [~,mdl_poles, b_mdl, ~]        = get_MP_components(H, 'll_mdl', params);
        [~,aic_poles, b_aic, ~]        = get_MP_components(H, 'll_aic', params);
        [~,map_poles, b_map, ~]        = get_MP_components(H, 'll_map', params);
        [~,evt_poles, b_evt, ~]        = get_MP_components(H, 'll_evt', params);

        SDD.Amplitudes.SV{ind,rep}    = b_SDD;
        SDD.Poles.SV{ind,rep}         = SDD_poles;
        SDD.NOF.SV(ind,rep)           = length(b_SDD);

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
        
        MAP.Amplitudes.LL{ind,rep}    = b_map;
        MAP.Poles.LL{ind,rep}         = map_poles;
        MAP.NOF.LL(ind,rep)           = length(b_map);

        EVT.Amplitudes.LL{ind,rep}    = b_evt;
        EVT.Poles.LL{ind,rep}         = evt_poles;
        EVT.NOF.LL(ind,rep)           = length(b_evt);
    end
    
    ind = ind + 1;

end

end

