function [M_opt_AIC, M_opt_MDL, M_opt_MAP, M_opt_EVT] = get_INFORMATION_CRITERIA_MP_model_order(H, M_max, params)

% number of parameter is pre-defined to be 4 - 
% initial phase, amplitude, damping factor and frequency
N_PARAMS = 4;

% params.true_poles = exp(-params.true_damps*params.dt + 1i*params.true_freqs*params.dt);
T  = params.interval_length;
dt = params.dt;
N = length(params.signal);

AIC = zeros(M_max,1);
MDL = zeros(M_max,1);
MAP = zeros(M_max,1);
EVT = zeros(M_max,1);

    for M = 1: M_max
        % Estimate the parameters for a given model order
        [est_poles, est_amps] = get_mp_est(H,M);
    
        % Compute the clean and estimated signals (assuming initial phases == 0)
        x_hat = sum( (est_amps)   .* (est_poles)   .^(0:T/dt),1) ;
        x     = sum( (params.true_amps.').* (params.true_poles.').^(0:T/dt),1) ;
    
        % Sanity test
        % x_emp = params.signal-params.noise;
        % norm(x- x_emp)
    
        % Compute log-likelihood
        logL = -N*log(pi*params.noise_variance) - (1/(params.noise_variance)) * norm(x_hat-x,2)^2;
    
        % Number of parameters
        n = N_PARAMS * M; 
        AIC(M) = n - logL;
        MDL(M) = n*log(N) - 2*logL;
        MAP(M) = 5*M*log(N)-2*logL;
        alpha = 0.05;
        EVT(M) = 2*M*(log(N) + log(log(N)) - log(3*(alpha^2)/pi)) -2*logL;
    end
    
    [~, M_opt_AIC] = min(AIC);
    [~, M_opt_MDL] = min(MDL);

    [~, M_opt_MAP] = min(MAP);
    [~, M_opt_EVT] = min(EVT);

end

function [est_poles, est_amps] = get_mp_est(H,M)

    X1 = H(:,1:end-1);
    X2 = H(:,2:end);

    [U, S, V] = svd(X1, 'econ');
    M = min(M, length(S));
    Ur = U(:, 1:M); 
    Sr = S(1:M, 1:M);
    Vr = V(:, 1:M);

    Ltilde =  Sr\(Ur')*X2*Vr;

    [W, D] = eig(Ltilde);
    est_poles = diag(D);

    left_mode  = Ur*Sr*W;
    right_mode = W \ Vr';       
    
    est_amps  = (left_mode(1,:).') .* right_mode(:,1); 

end