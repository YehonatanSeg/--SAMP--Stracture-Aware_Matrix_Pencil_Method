function [params, TIMES, EST] = gen_data_amps_comparisons_MP_SAMP(noise_distribution, Nreps, damps)
%% General settings
params.RANK_TYPE = 'numeric';
params.NOSIE_TYPE = noise_distribution;
params.Nrep      = Nreps;        
params.SNR_range = 10;
params.interval_length_range = 5:5:35;

params.dt = 1/10;
params.interval_length = 10;
N = 1 +params.interval_length / params.dt;


%% Signals settings
w1 = 2;
params.true_freqs  = [w1 + 1.5*pi/(params.dt*N), w1 ] ; 
params.true_damps  = damps; 
params.true_amps   = [1, 1]; 
params.true_phases = [0, 0];

% Defining the data structures
new_amps_times  = zeros(length(params.interval_length_range),params.Nrep);
MP_amps_times   = zeros(length(params.interval_length_range),params.Nrep);

new_amps_est1  = zeros(length(params.interval_length_range),params.Nrep);
new_amps_est2  = zeros(length(params.interval_length_range),params.Nrep);

MP_amps_est1   = zeros(length(params.interval_length_range),params.Nrep);
MP_amps_est2   = zeros(length(params.interval_length_range),params.Nrep);

% Main calculation
SNR = params.SNR_range;
k=1;
for intLength = params.interval_length_range
    for rep=1:params.Nrep
        disp(['[signal length, rep] = [', num2str(intLength),', ',num2str(rep),']' ])
        [params.signal, params.noise, params.noise_variance]  = generate_sum_of_k_complex_exp_with_add_noise(params.true_freqs, params.true_damps, params.true_phases, params.true_amps, params.dt, intLength, SNR, params.NOSIE_TYPE);

        params.pencil_parameter = round(length(params.signal) /3);

        Y = create_hankel_matrix(params.signal, params.pencil_parameter);
        Y0 = Y(:,1:end-1);
        Y1 = Y(:,2:end);

        [U, S, V] = svd(Y0, 'econ');
        r = rank(Y0); 
        Ur = U(:, 1:r); 
        Sr = S(1:r, 1:r);
        Vr = V(:, 1:r);
    
        Ltilde =  Sr\(Ur')*Y1*Vr;

        [Q, D] = eig(Ltilde);
        lambda = diag(D);
    
        tic;
        left_mode  = Ur*Sr*Q;
        right_mode = Q \ Vr';  
        b_new  = (left_mode(1,:).') .* right_mode(:,1);  
        t_new = toc;

        % the classical MP way 
        tic
        Vander = create_vander(lambda, sum(size(Y))-1).';
        y = (params.signal).';
        b_MP = pinv(Vander)*y;
        t_MP = toc;

        TIMES.SAMP(k,rep)  = t_new;
        TIMES.MP(k,rep)    = t_MP;

        [EST.SAMP_est1(k,rep), EST.SAMP_est2(k,rep)] = get_estimation(params.true_amps, abs(b_new));
        [EST.MP_est1(k,rep), EST.MP_est2(k,rep)]     = get_estimation(params.true_amps, abs(b_MP));

    end
    k=k+1; 
end
end

%% Aux Function

% true_amps are assumed to be sorted in descending order
% and number of estimated freqs is assumed to be greater than 2 !
% this is important when using the sort function
function [est1, est2] =  get_estimation(true_amps, estim_freqs)

    [est_1, est_2] = closestNumbers(true_amps(1), true_amps(2), estim_freqs);
    
    % If only one frequency has been estimated
    if ~sum(est_2) 
        est_2 = inf;
    elseif ~sum(est_1)
        est_1 = inf;
    end

    sorted_estim_amps = sort([est_1, est_2], 'descend');

    est1 = sorted_estim_amps(1); 
    est2 = sorted_estim_amps(2); 

    
end

