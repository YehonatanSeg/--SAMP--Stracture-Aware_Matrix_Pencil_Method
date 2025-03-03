function epsilon = compute_SAMP_signal_detection_features(left_mode, poles, params)

% Extract the relevant information from the left mode
Nfft = 2^(nextpow2(size(left_mode,1)));
y = fft(left_mode-mean(left_mode,1), Nfft, 1);
y = abs(y); 
[max_mode_peaks, ~]  = max(y);
epsilon = (max_mode_peaks') ;  % Signal detection features by Eq (41)
freqs = angle(poles)/params.dt;
diff_matrix = abs(freqs - freqs.').^2;
freqs_sum_diff = sum(diff_matrix, 2);
epsilon = epsilon./freqs_sum_diff;  % Normalized signal detection features by Eq (42)

end