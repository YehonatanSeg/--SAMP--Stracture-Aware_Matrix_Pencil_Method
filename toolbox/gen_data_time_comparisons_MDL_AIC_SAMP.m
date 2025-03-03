function [Poles, Times,NOF, params] = gen_data_time_comparisons_MDL_AIC_SAMP(noise_dist, Nreps, damps)
%% General settings

%% Settings
params.RANK_TYPE = 'effective';
params.NOSIE_TYPE = noise_dist;
params.Nrep = Nreps;
params.SNR_RANGE = 10;
params.FREQ_DIST = 1.5;
params.INTERVAL_LENGTH_RANGE = 5:5:35;
SNR = params.SNR_RANGE;
params.dt = 1/10;
params.true_amps   = [1, 1];
params.true_phases = [0, 0];
w1 = 2; 
 
ind = 1;
for intLength = params.INTERVAL_LENGTH_RANGE 
    
    params.interval_length = intLength; 
    N = 1 +params.interval_length / params.dt;
    params.true_freqs  = [w1 + params.FREQ_DIST*pi/(params.dt*N), w1 ] ; 
    params.true_damps  = damps;
    params.true_poles = exp(-params.true_damps*params.dt + 1i*params.true_freqs*params.dt);
    
    for rep = 1:params.Nrep
        disp(['Dampings = [', num2str(damps), ']', ' interval length = ', num2str(intLength), ' SNR = ', num2str(SNR), ' rep = ', num2str(rep)])

        [params.signal, params.noise, params.noise_variance]  = generate_sum_of_k_complex_exp_with_add_noise(params.true_freqs, params.true_damps, params.true_phases, params.true_amps, params.dt, intLength, SNR, params.NOSIE_TYPE);
        params.pencil_parameter = round(length(params.signal) /3);
        H = create_hankel_matrix(params.signal,params.pencil_parameter);

        
        % SAMP (our)
        tic;
        [left_mode, mp_eigenvalues, b_MP, b_new] = get_MP_components(H, params.RANK_TYPE, []);
        epsilon = compute_SAMP_signal_detection_features(left_mode, mp_eigenvalues, params);
        t_eff_MP = toc;

        tic;
        ind_SAMP_GAP    = get_true_idx(epsilon, 'GAP');
        t_SAMP_GAP = toc;
        Times.SAMP_gap(ind, rep)  = t_eff_MP+t_SAMP_GAP;

        tic;
        ind_SAMP_KMEANS = get_true_idx(epsilon, 'kmeans');
        t_SAMP_KMEANS = toc;
        Times.SAMP_kmeans(ind, rep)  = t_eff_MP+t_SAMP_KMEANS;

        Poles.SAMP_gap{ind, rep}     = get_poles_est(mp_eigenvalues(ind_SAMP_GAP), params.true_poles);
        Poles.SAMP_kmeans{ind, rep}  = get_poles_est(mp_eigenvalues(ind_SAMP_KMEANS), params.true_poles);
        
        NOF.SAMP_gap(ind, rep)    = length(ind_SAMP_GAP);
        NOF.SAMP_kmeans(ind, rep) = length(ind_SAMP_KMEANS);

    
        % AIC \ MDL
        tic;
        [~,mdl_poles, b_mdl, ~]        = get_MP_components(H, 'll_mdl', params);
        t_MDL = toc;

        tic;
        [~,aic_poles, b_aic, ~]        = get_MP_components(H, 'll_aic', params);
        t_AIC = toc;

        Times.MDL(ind, rep)    = t_MDL;
        Times.AIC(ind, rep)    = t_AIC;
        
        Poles.MDL{ind,rep} = get_poles_est(mdl_poles, params.true_poles);
        Poles.AIC{ind,rep} = get_poles_est(aic_poles, params.true_poles);
        
        NOF.MDL(ind, rep) = length(b_mdl);
        NOF.AIC(ind, rep) = length(b_aic);

    end
    
    ind = ind + 1;

end

end

%% AUX functions 
function est_poles =  get_poles_est(poles, true_poles)
    
        [est1, est2] = closestNumbers(true_poles(1), true_poles(2), poles);
        
        % If only one frequency has been estimated
        if isempty(est2) 
            est2 = inf;
        elseif isempty(est1)
            est1 = inf;
        end
        
        % sort by freqs in descending order
        if angle(est1) > angle(est2)
            est_poles(1) = est1;
            est_poles(2) = est2;
        else
            est_poles(1) = est2;
            est_poles(2) = est1;
        end

    
end
   
 
