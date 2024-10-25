% y: the measured signal (noisy signal)
% b1, b2: amplitudes  
% phi1, phi2: the initial phases
% alpha1, alpha2: damping factors
% theta1, theta2: frequencies 
% var: variance of the noise
% N: number of samples
% dt: sampling rate
function CRB = CRB_unbiased_full_model(b1, b2, phi1, phi2, alpha1, alpha2, theta1, theta2, var, N, dt)
  
    % Initialize the Fisher Information (FIM) Matrix 
    FIM = zeros(8, 8); 
    
    s1 = (-alpha1 + 1j * theta1);
    s2 = (-alpha2 + 1j * theta2);

    % Compute the FIM components
    for n = 0:N-1

        ndt = n * dt; 

        c_b1     = -exp(1j * phi1) * exp(s1 * ndt) ; 
        c_b2     = -exp(1j * phi2) * exp(s2 * ndt); 

        c_phi1   = -1j * b1 * exp(1j * phi1) * exp(s1 *ndt) ; 
        c_phi2   = -1j * b2 * exp(1j * phi2) * exp(s2 * ndt) ; 

        c_alpha1 = ndt * b1 * exp(1j * phi1) * exp(s1 * ndt) ; 
        c_alpha2 = ndt * b2 * exp(1j * phi2) * exp(s2 * ndt) ; 

        c_theta1 = -1j * ndt * b1 * exp(1j * phi1) * exp(s1 * ndt) ; 
        c_theta2 = -1j * ndt * b2 * exp(1j * phi2) * exp(s2 * ndt) ; 

        % Gradient vector
        gradient = [c_b1; c_b2; c_phi1; c_phi2; c_alpha1; c_alpha2; c_theta1; c_theta2];
        
        % Update the FIM
        FIM = FIM + real(gradient * gradient') ;
    end

    % Scaling
    FIM =  FIM * (2/var); 

%%
    
    % Compute the inverse of the Fisher Information Matrix
    inv_FIM = inv(FIM);

    % The CRB for the params are the diagonal elements of the inverse FIM 
    CRB.b1     = inv_FIM(1, 1);
    CRB.b2     = inv_FIM(2, 2);    
    CRB.phi1   = inv_FIM(3, 3);
    CRB.phi2   = inv_FIM(4, 4);    
    CRB.alpha1 = inv_FIM(5, 5);
    CRB.alpha2 = inv_FIM(6, 6);
    CRB.theta1 = inv_FIM(7, 7);
    CRB.theta2 = inv_FIM(8, 8);
end

