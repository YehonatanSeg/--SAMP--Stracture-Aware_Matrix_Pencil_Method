function epsilon = compute_SAMP_plusplus_signal_detection_features(left_mode, poles, params)

% Extract the relevant information from the left mode

% Mode energy term by Eq (18) 
Nfft = 2^(nextpow2(size(left_mode,1)));
y = fft(left_mode-mean(left_mode,1), Nfft, 1);
y = abs(y); 
[max_mode_peaks, ~]  = max(y);
rho = (max_mode_peaks') ; 
freqs = angle(poles)/params.dt;
diff_matrix = abs(freqs - freqs.').^2;
freqs_sum_diff = sum(diff_matrix, 2);
rho = rho./freqs_sum_diff;  


% Mode correlations term by Eq (19)
cov_left_modes = abs(cov(left_mode));
off_diag_cov_left_modes = cov_left_modes;
off_diag_cov_left_modes(logical(eye(size(off_diag_cov_left_modes)))) = NaN;
sigma = std(off_diag_cov_left_modes,[],2,"omitmissing");
    
    
% Signal detection features by Eq (17) 
alpha = 0.5;
epsilon = alpha* (rho/max(rho)) + (1-alpha)*(sigma/max(sigma)); 


end