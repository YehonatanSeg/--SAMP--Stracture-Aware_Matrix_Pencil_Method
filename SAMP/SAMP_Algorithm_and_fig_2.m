addpath(fullfile(fileparts(pwd), 'toolbox'));
clear;
close all;
%% Simulating a sum of complex exponentials and additive noise

% Settings for simulating sum of complex exponentials by Eq (2)
params.dt = 1/10;
params.interval_length = 7;
N = 1 +params.interval_length / params.dt;

w1 = 2; 

params.true_freqs  = [w1 + 2*pi/(params.dt*N), w1] ; 
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
[left_modes, right_modes, mp_poles] = get_MP_components(Y); 

% Parameter Estimation by step 2:
mp_freqs  = imag(log(mp_poles))/params.dt;
mp_damps  = real(log(mp_poles)) /params.dt;
b_tilde   = (left_modes(1,:).') .* right_modes(:,1); % Efficient amplitude estimation by Eq (44)

% Model Order Detection by step 3:

%  Normalized similaruty measure by Eq (36)
P = compute_similarity_measure(left_modes, params);

% Signal detection features by Eq (40)
epsilon = max(P  , [] , 2 ); 

% Normalized signal detection features by Eq (43)
diff_matrix = (mp_poles ./ mp_poles.'); 
damps_sum_diff = sum(abs(diff_matrix).^2 , 1).^2;
epsilon = epsilon ./ damps_sum_diff(:);
epsilon = epsilon/max(epsilon);


% Feature binary partition into two distinct subsets by
% defining the signal-related subset as the subset of
% features satisfying Eq. (46):
bt_i        =  abs(b_tilde);              
N_L         = size(left_modes,1);
A           = @(z) (z .^ (0:N_L-1)).';
c           = 10*sqrt(N_L)./ vecnorm(A(mp_poles)); 
ISR         = c(:) ./ bt_i ;
P_low_bound = ( (1-ISR) ./ (1+ISR) ).^2 ;
signal_ind  = find( epsilon(:) > P_low_bound(:) );

% model order estimation
M_tilde     = length(signal_ind); 

% Parameter Selection by step 4:
est_freqs   = mp_freqs(signal_ind);
est_damps   = mp_damps(signal_ind);
est_amps    = b_tilde(signal_ind);

%% Figs
%Fig. (2) in the paper

fig = figure();  hold on;

all_dist = abs(mp_freqs - params.true_freqs).^2;
pairs = matchpairs(all_dist,1e16);
xline(pairs(:,1), 'LineWidth',3,'color','red'); 
plot(epsilon,'o', 'LineWidth',2); 
plot(P_low_bound, 'x', 'color','black', 'LineWidth',2);
xlabel('Mode index','FontSize',14)
ylabel('$P_i(z_i^*)$', 'Interpreter','latex','FontSize',14)

%% ======================= AUX FUNCTIONS ======================= %%
function P = compute_similarity_measure(left_modes, params)
    % Define the Grid of frequencies x dampings
    N_grid = 2^(nextpow2(size(left_modes,1)));
    thetaGrid = linspace(-pi, pi, N_grid);
    damp_grid = linspace(0, 1, N_grid)';
    lambda_grid = exp(-(damp_grid*params.dt) + 1i*thetaGrid*params.dt);
    
    left_modes_n = left_modes./vecnorm(left_modes);
    N = size(left_modes_n,1);
    
    % --- build Vandermonde for *all* grid points ----------------------------
    zvec   = lambda_grid(:).';                     
    Agrid  =  zvec.^((0:N-1).') ;
    denom  = sum(abs(Agrid).^2,1);      
    
    % --- numerator for every mode & grid point in one shot ------------------
    num    = abs( left_modes_n' * Agrid ).^2;        
    
    % --- Compute the normalized similarity measure --------------------------
    P = num ./ denom;
end

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