%% Simulating a sum of complex exponentials and additive noise

% Settings for simulating sum of complex exponentials by Eq (2)
params.dt = 1/10;
params.interval_length = 10;
N = 1 +params.interval_length / params.dt;
w1 = 2;
params.true_freqs  = [w1 + 1.5*pi/(params.dt*N), w1 ] ; 
params.true_damps  = [0.0, 0.0];
params.true_amps   = [1, 1];
params.true_phases = [0, 0];
params.true_poles = exp(-params.true_damps*params.dt + 1i*params.true_freqs*params.dt);
SNR = 10;

% Generating the noisy signal y(n) by Eq (1)
[params.signal, params.noise, params.noise_variance]  = generate_sum_of_k_complex_exp_with_add_noise(params.true_freqs, params.true_damps, params.true_phases, params.true_amps, params.dt, params.interval_length, SNR, 'normal');

% Defining the Pencil parameter ( L = round(N/3) )
params.pencil_parameter = round(length(params.signal)/3);

% Creating the noisy Hankel matrix by Eq (3)
Y = create_hankel_matrix(params.signal,params.pencil_parameter);
%% The SAMP Algorithm (Algorithm 1)

%Modes and Eigenvalues Computation by step 1:
[left_mode, right_mode, lambda_mp] = get_MP_components(Y, 'effective',[]); % ll_mdl, ll_aic

% Parameter Estimation by step 2:
mp_freqs  = imag(log(lambda_mp))/params.dt;
mp_damps  = real(log(lambda_mp)) /params.dt;
b_tilde   = (left_mode(1,:).') .* right_mode(:,1); % Efficient amplitude estimation by Eq (43)

% Model Order Detection by step 3:
Nfft = 2^(nextpow2(size(left_mode,1)));
y = fft(left_mode-mean(left_mode,1), Nfft, 1);
y = abs(y); 
[max_mode_peaks, ~]  = max(y);
epsilon = (max_mode_peaks') ;  % Signal detection features by Eq (41)
freqs = angle(lambda_mp)/params.dt;
diff_matrix = abs(freqs - freqs.').^2;
freqs_sum_diff = sum(diff_matrix, 2);
epsilon = epsilon./freqs_sum_diff;  % Normalized signal detection features by Eq (42)

% Deviding the signal detection features into two distinct groups using
% k-means algorithm and selecting the signal related indices 
idx            = kmeans(epsilon,2,'Replicates', 1);
avg(1)         = mean(epsilon( idx == 1 ));
avg(2)         = mean(epsilon( idx == 2 ));
[~,true_group] = max(avg);
signal_ind     = find(idx == true_group);
M_tilde        = length(signal_ind); % model order estimation

% Parameter Selection by step 4:
est_freqs = mp_freqs(signal_ind);
est_damps = mp_damps(signal_ind);
est_amps  = b_tilde(signal_ind);

%% Figs
%Fig. (2) in the paper
plot_FFT_MSMP(left_mode, params, signal_ind)











%% ======================= AUX FUNCTIONS ======================= %%

function [left_mode, right_mode, lambda] = get_MP_components(Y, rank_type, params)

        Y1 = Y(:,1:end-1);
        Y2 = Y(:,2:end);

        [U, S, V] = svd(Y1, 'econ');    
        r = get_rank(Y1, rank_type, params);
        
        Ur = U(:, 1:r); 
        Sr = S(1:r, 1:r);
        Vr = V(:, 1:r);
    
        Ltilde =  Sr\(Ur')*Y2*Vr;

        [W, D] = eig(Ltilde);
        lambda = diag(D);
    
        left_mode  = Ur*Sr*W;
        right_mode = W \ Vr';       
        
end

function r = get_rank(A, rank_type, params)
    
    [~, S, ~] = svd(A,'econ');
    
    if strcmp(rank_type, 'numeric')
        r = rank(A);
    elseif strcmp(rank_type, 'effective')
        r = round(erank(A));
    elseif strcmp(rank_type, 'significant')
        [r, ~] = get_rank_by_significant_digits(diag(S), 2); 
    elseif strcmp(rank_type, 'gap')
        [r, ~] = get_rank_by_singular_value_gap(diag(S));
    elseif strcmp(rank_type, 'kmeans')
        [r, ~] = get_rank_by_singular_value_kmeans(diag(S));
    elseif strcmp(rank_type, 'll_mdl')    
        [~, r] = get_AIC_MDL_MP_model_order(A, params);
    elseif strcmp(rank_type, 'll_aic')    
        [r, ~] = get_AIC_MDL_MP_model_order(A, params);
    end


end
