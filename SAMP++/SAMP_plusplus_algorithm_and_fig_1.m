addpath(fullfile(fileparts(pwd), 'toolbox'));
clear;
close all;
%% Simulating a sum of complex exponentials and additive noise

% Settings for simulating sum of complex exponentials by Eq (1)
params.dt = 1/15;
params.interval_length = 10;
N = 1 +params.interval_length / params.dt;
w1 = 2;
w2 = 4;
params.true_freqs  = [w2 + 1.5*pi/(params.dt*N), w2 ,w1 + 1.5*pi/(params.dt*N), w1] ; 
params.true_damps  = [0.0, 0.0,0.0, 0.0];
params.true_amps   = [1, 1, 1, 1];
params.true_phases = [0, 0,0, 0];

params.true_poles = exp(-params.true_damps*params.dt + 1i*params.true_freqs*params.dt);
SNR = 20;

% Generating the noisy signal y(n) by Eq (1)
[params.signal, params.noise, params.noise_variance]  = generate_sum_of_k_complex_exp_with_add_noise(params.true_freqs, params.true_damps, params.true_phases, params.true_amps, params.dt, params.interval_length, SNR, 'normal');

% Defining the Pencil parameter ( L = round(N/3) )
params.pencil_parameter = round(length(params.signal)/3);

% Creating the noisy Hankel matrix by Eq (2)
Y = create_hankel_matrix(params.signal,params.pencil_parameter);

%% The SAMP++ Algorithm (Algorithm 1)

% Modes and Eigenvalues Computation by step 1:
[left_mode, right_mode, lambda_mp] = get_MP_components(Y); % ll_mdl, ll_aic

% Parameter Estimation by step 2:
mp_freqs  = imag(log(lambda_mp))/params.dt;
mp_damps  = real(log(lambda_mp)) /params.dt;
b_tilde   = (left_mode(1,:).') .* right_mode(:,1); % Efficient amplitude estimation by Eq (43) in SAMP paper

% Model Order Detection by step 3:

% Mode energy term by Eq (18) 
Nfft = 2^(nextpow2(size(left_mode,1)));
y = fft(left_mode-mean(left_mode,1), Nfft, 1);
y = abs(y); 
[max_mode_peaks, ~]  = max(y);
rho = (max_mode_peaks') ; 
freqs = angle(lambda_mp)/params.dt;
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


% Deviding the signal detection features into two distinct groups using
% k-means algorithm and selecting the signal related indices 
idx            = kmeans(epsilon,2,'Replicates', 5);
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
%Fig. 1 in SAMP++
plot_modes_cov(left_mode, signal_ind)




%% ======================= AUX FUNCTIONS ======================= %%

function [left_mode, right_mode, lambda] = get_MP_components(Y)

        Y0 = Y(:,1:end-1);
        Y1 = Y(:,2:end);

        % SVD of Y0 by Eq (4)
        [U, S, V] = svd(Y0, 'econ'); 

        % Compute the effective rank by Eq (44)
        S_norm = sum(abs(diag(S)));
        p = diag(S)/S_norm;
        H = -sum(p.*log(p));
        r = round(exp(H));

        % Truncating the SVD of Y0
        Ur = U(:, 1:r); 
        Sr = S(1:r, 1:r);
        Vr = V(:, 1:r);
        
        % Computing the matrix A by Eq (6)
        A =  Sr\(Ur')*Y1*Vr;
        
        % Eigendecomposition of A by Eq (7)
        [Q, D] = eig(A);
        lambda = diag(D);
    
        left_mode  = Ur*Sr*Q;
        right_mode = Q \ Vr';       
        
end

function plot_modes_cov(left_mode, true_ind)


% Computing the absolute value of the sample covariance matrix
centered_left_modes = left_mode - mean(left_mode);
cov_left_modes = abs(centered_left_modes'*centered_left_modes);

% Keep only the off diagonal elements of the sample covariance matrix
off_diag_cov_left_modes = cov_left_modes;
off_diag_cov_left_modes(logical(eye(size(off_diag_cov_left_modes)))) = NaN;
validIdx = ~isnan(off_diag_cov_left_modes); 

% Randomly select the noise related indices
noise_ind = 1:length(cov_left_modes);
noise_ind(true_ind) = [];
two_randomly_noise_ind = noise_ind(randperm(length(noise_ind), 2));


% Select two signal-related row of the cov matrix 
C_true_1 = off_diag_cov_left_modes(true_ind(1),validIdx(true_ind(1),:));
C_true_2 = off_diag_cov_left_modes(true_ind(2),validIdx(true_ind(2),:));

% Select two noise-related row of the cov matrix 
C_noise_1 = off_diag_cov_left_modes(two_randomly_noise_ind(1),validIdx(two_randomly_noise_ind(1),:));
C_noise_2 = off_diag_cov_left_modes(two_randomly_noise_ind(2),validIdx(two_randomly_noise_ind(2),:));

% Plot Fig 1
fig = figure();
fig.WindowState = 'maximized'; 
pause(1)


% These setting are for placing 4 subfigures in one column..
set(0,'DefaultTextFontName','Times',...
'DefaultTextFontSize',75,...
'DefaultAxesFontName','Times',...
'DefaultAxesFontSize',75,...
'DefaultLineLineWidth',4.5,...
'DefaultLineMarkerSize',15);

% Set max ylim, xlim for all figures
ylim_max = max(max(max(C_true_1), max(C_true_2)), max(max(C_noise_1), max(C_noise_2)));
x_lim_max = size(left_mode,2);

% Signal row 1
subplot(2,2,1);
plot(find(validIdx(true_ind(1),:)), C_true_1); 
ylim([0,ylim_max]);
xlim([0,x_lim_max]);
ylabel('$|[\mathbf{C}]_{i,k}|$','Interpreter','latex')
xlabel('Indices')

% Signal row 2
subplot(2,2,2);
plot(find(validIdx(true_ind(2),:)), C_true_2); 
ylim([0,ylim_max]);
xlim([0,x_lim_max]);
ylabel('$|[{\mathbf{C}}]_{i,k}|$','Interpreter','latex')
xlabel('Indices')

% Noise row 1
subplot(2,2,3);
plot(find(validIdx(two_randomly_noise_ind(1),:)), C_noise_1)
ylim([0,ylim_max]);
xlim([0,x_lim_max]);
ylabel('$|[{\mathbf{C}}]_{i,k}|$','Interpreter','latex')
xlabel('Indices')

% Noise row 2
subplot(2,2,4);
plot(find(validIdx(two_randomly_noise_ind(2),:)), C_noise_2)
ylim([0,ylim_max]);
xlim([0,x_lim_max]);
ylabel('$|[{\mathbf{C}}]_{i,k}|$','Interpreter','latex')
xlabel('Indices')

end

