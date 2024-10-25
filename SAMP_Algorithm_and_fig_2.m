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
[left_mode, right_mode, lambda_mp] = get_MP_components(Y); % ll_mdl, ll_aic

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
plot_modes_FFT(left_mode, params, signal_ind)




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


function plot_modes_FFT(left_mode, params, true_ind)

Nfft = 2^(3+nextpow2(size(left_mode,1)));
y = fft(left_mode-mean(left_mode,1),Nfft, 1);
y = fft(left_mode,Nfft, 1);

y = abs(y);
global_max = max(y,[],'all');
y = y/global_max;
frequencies=(0:Nfft-1)*(2*pi/(params.dt*Nfft));

% re-ordering the indices according to the true-frequencies
[~, freqs_ind] = max(y(:,true_ind));
est_freqs = frequencies(freqs_ind);
dist = abs(est_freqs' - params.true_freqs);
if sum(diag(dist)) > sum(diag(flip(dist)))
    true_ind_ord(1) = true_ind(2);
    true_ind_ord(2) = true_ind(1);
    true_ind = true_ind_ord;
end

% These setting are for placing 4 subfigures in one column..
set(0,'DefaultTextFontName','Times',...
'DefaultTextFontSize',75,...
'DefaultAxesFontName','Times',...
'DefaultAxesFontSize',75,...
'DefaultLineLineWidth',4.5,...
'DefaultLineMarkerSize',15);

fig1 = figure(); hold on;
fig1.WindowState = 'maximized';  
plot(frequencies,y(:,true_ind(1))); 
ylim([0,1]);
ylabel('Coefficients')
xlabel('$\theta$', 'Interpreter', 'latex')
plot([params.true_freqs(1) params.true_freqs(1)],[0 1], '-.','color', 'black');
legend('$|\mathcal{F}(\widetilde{\mathbf{Z}}_L^{1})|$', '$\theta_1$', 'Interpreter', 'latex')

path = 'C:\Users\yehon\Desktop\Phd Studies\Projects\SAMP\Code\Images';
name = 'DFT_signal_modes_1.png';
exportgraphics(fig1,[path,'\', name],'Resolution',300) 

fig2 = figure(); hold on;
fig2.WindowState = 'maximized';
plot(frequencies,y(:,true_ind(2))); 
ylim([0,1]);
ylabel('Coefficients')
xlabel('$\theta$', 'Interpreter', 'latex')
plot([params.true_freqs(2) params.true_freqs(2)],[0 1], '-.','color', 'black');
legend('$|\mathcal{F}(\widetilde{\mathbf{Z}}_L^{2})|$', '$\theta_2$', 'Interpreter', 'latex')

name = 'DFT_signal_modes_2.png';
exportgraphics(fig2,[path,'\', name],'Resolution',300) 

fig3 = figure();
fig3.WindowState = 'maximized';
plot(frequencies,y(:,max(true_ind)+1))
ylim([0,1]);
ylabel('Coefficients')
xlabel('$\theta$', 'Interpreter', 'latex')
legend('$|\mathcal{F}(\widetilde{\mathbf{Z}}_L^{3})|$', 'Interpreter', 'latex')

name = 'DFT_noise_modes_1.png';
exportgraphics(fig3,[path,'\', name],'Resolution',300) 

fig4 = figure();
fig4.WindowState = 'maximized';
plot(frequencies,y(:,20))
ylim([0,1]);
ylabel('Coefficients')
xlabel('$\theta$', 'Interpreter', 'latex')
legend('$|\mathcal{F}(\widetilde{\mathbf{Z}}_L^{20})|$', 'Interpreter', 'latex')

name = 'DFT_noise_modes_2.png';
exportgraphics(fig4,[path,'\', name],'Resolution',300) 

end