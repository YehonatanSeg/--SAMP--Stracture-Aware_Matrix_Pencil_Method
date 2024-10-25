function [f, noise, noise_var] = generate_sum_of_k_complex_exp_with_add_noise(freqs, alpha, phases, amps, dt, T, SNR, noise_type)

if length(freqs) ~= length(alpha)
    error('Number of freqs and dumps is different')
end

if length(freqs) ~= length(amps)
    error('Number of freqs and amps is different')
end

if any(freqs >= 1/(2*dt))
    error('Sampling rate is too low w.r.t Nyquist')
end

if freqs ~= sort(freqs, 'descend')
        error('Freqs should be in decending order')
end

if amps ~= sort(amps, 'descend')
        error('Amps should be in decending order')
end

t = 0:dt:T;

s = zeros(length(freqs));
for i=1: length(freqs)
    s(i)   = -alpha(i)  + 1i*freqs(i);
    h(i,:) = amps(i) * exp(s(i) .* t + 1i*phases(i)) ;
end



% fc = sum(amps.*exp((-alpha + 1i*freqs)*t + 1i*phases),1);
fc = sum(h,1);

noise = [];
noise_var=[];
f = fc;
if ~isempty(SNR)
    
    % white, complex Gaussian noise 
    if strcmp(noise_type, 'normal')
        n = randn(size(fc)) + 1i*randn(size(fc));
    elseif strcmp(noise_type, 'uniform')
        n = rand(size(fc))  + 1i*rand(size(fc));
    elseif strcmp(noise_type, 'bi_normal')
        n = randsn(size(fc), 0, 2, 1, 0.1) + 1i*randsn(size(fc), 0, 2, 1, 0.1);
    elseif strcmp(noise_type, 't_student')
        n = trnd(1, size(fc)) + 1i*trnd(1, size(fc));
    end

    Es = sum(amps.^2);
    En = var(n);

    SNR_linear = 10^(SNR/ 10);

    alpha = sqrt( Es/ ( En * SNR_linear ) );
    f     = fc+alpha*n;
    noise = f-fc;
   

    noise_var = var(noise);
end


end

function samples = randsn(N, mu1, mu2, sigma1, sigma2)
% Parameters for mixture of normals
% mu1 = 0; sigma1 = 1; % for the main normal distribution
% mu2 = 3; sigma2 = 2; % for the tail distribution
p = 0.85; % probability for the main distribution

% Generate random numbers
r = rand(N);
samples = (r < p).*normrnd(mu1, sigma1, N) + ...
          (r >= p).*normrnd(mu2, sigma2, N);

end