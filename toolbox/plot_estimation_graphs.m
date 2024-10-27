function plot_estimation_graphs(params, CRB, AIC, MDL, GAP,EFF, KMEANS, x_axis_param)
%% GENERAL CONFIG
METHODS = {'AIC', 'MDL', 'GAP', 'EFF', 'KMEANS'}; 
SNR_RANGE_GRAPHS = [0,30];
set_IEEE_TSP_default_figure_settings()
X_AXIS.x_axis_param = x_axis_param;

if strcmp(X_AXIS.x_axis_param, 'SNR')
    X_AXIS.x = params.SNR_RANGE;
    X_AXIS.x_axis_range = SNR_RANGE_GRAPHS;
    X_AXIS.x_axis_label = 'SNR [dB]';
elseif strcmp(X_AXIS.x_axis_param, 'FREQS_DIST')
    X_AXIS.x = params.FREQ_DIST_RANGE;
    X_AXIS.x_axis_range = [params.FREQ_DIST_RANGE(1), params.FREQ_DIST_RANGE(end)];
    X_AXIS.x_axis_label = '$|\Delta\theta|$ [rad/s]';
elseif strcmp(X_AXIS.x_axis_param, 'SAMPLES')
    X_AXIS.x = params.INTERVAL_LENGTH_RANGE / params.dt;
    X_AXIS.x_axis_range = [X_AXIS.x(1), X_AXIS.x(end)];
    X_AXIS.x_axis_label = 'N';
end

%% Best single param (frequency, damping) estimation graphs

% Averaging over all Repetitions
%  As we have inf when estimating only a single pole we cannot use RMSE or
%  any mean, this we use median.
Single_POLE_RMSE = [];
Single_DAMP_RMSE = [];
Single_FREQ_RMSE = [];
% producing median resutls for Figure 2 and 3 (undapmed and damped figures
% for the combination of two poles)

% Define the desired component
comp_num = 1;

for k=1:length(METHODS)
    method = eval(METHODS{k});
    fields = fieldnames(method.Poles);
    for i =1:length(fields)
        
        % RMSE COMPUTATIONS
        [curr_pole_error, freq_est, curr_frequency_error, ~, curr_damping_error]  = get_best_params_error(method.Poles.(fields{i}), params, comp_num, 'rmse', X_AXIS.x_axis_param) ;
        Single_POLE_RMSE.(METHODS{k}).(fields{i})      = sqrt(mean(curr_pole_error,2));
        Single_FREQ_RMSE.(METHODS{k}).(fields{i})      = sqrt(mean(curr_frequency_error,2));
        Single_DAMP_RMSE.(METHODS{k}).(fields{i})      = sqrt(mean(curr_damping_error,2));

        % VARIANCE COMPUTATIONS
        [curr_pole_error, ~, curr_frequency_error, ~, curr_damping_error]  = get_best_params_error(method.Poles.(fields{i}), params, comp_num, 'var', X_AXIS.x_axis_param) ;
        Single_POLE_VAR.(METHODS{k}).(fields{i})       = var(curr_pole_error,[],2);
        Single_FREQ_VAR.(METHODS{k}).(fields{i})       = var(curr_frequency_error,[],2);
        Single_DAMP_VAR.(METHODS{k}).(fields{i})       = var(curr_damping_error,[],2);

        % BIAS COMPUTATIONS
        [curr_pole_error, ~, curr_frequency_error, ~, curr_damping_error]  = get_best_params_error(method.Poles.(fields{i}), params, comp_num, 'bias', X_AXIS.x_axis_param) ;
        Single_POLE_BIAS.(METHODS{k}).(fields{i})       = abs(mean(curr_pole_error,2));
        Single_FREQ_BIAS.(METHODS{k}).(fields{i})       = abs(mean(curr_frequency_error,2));
        Single_DAMP_BIAS.(METHODS{k}).(fields{i})       = abs(mean(curr_damping_error,2));

    end
end

%% FREQUENCY PLOTS

% Produce RMSE plots
freqs_RMSE_figs = produce_single_param_plots(Single_FREQ_RMSE, [], CRB.theta1, params, comp_num, 'RMSE', 'freqs', X_AXIS);

% Produce bais plots
freqs_bias_figs = produce_single_param_plots(Single_FREQ_BIAS, [], [], params, comp_num, 'bias', 'freqs', X_AXIS);

%% parameter distribution graphs
if strcmp(X_AXIS.x_axis_param, 'SNR')
% Define the desired component
SNR_level = 10;
ind = (X_AXIS.x == SNR_level);

ALL_EST_FREQS = [];
for k=1:length(METHODS)
    method = eval(METHODS{k});
    fields = fieldnames(method.Poles);
    for i =1:length(fields)    

        % freqs, damping and amplidues of all reps at a specific SNR/N/d_theta level 
        all_poles = method.Poles.(fields{i})(ind,:);
        all_poles =   vertcat(all_poles{:});
        ALL_EST_FREQS.(METHODS{k}).(fields{i})  = imag(log(all_poles)) / params.dt;
    end
end
 
%% Plots of the estimated parameters distribution

set_IEEE_TSP_default_figure_settings()
BW_freqs = abs(mean(diff(params.true_freqs)))/50;

% subplot(2,2,1); hold on
freqs_fig_MDL = figure(); hold on; 
freqs_fig_MDL.WindowState = 'maximized';
plot_histo(ALL_EST_FREQS.MDL.LL, [0.4660 0.6740 0.1880], 'MDL', ALL_EST_FREQS.KMEANS.SAMP, [0 0.4470 0.7410], 'SAMP (ours)' , BW_freqs, params.true_freqs)

% subplot(2,2,2); hold on
freqs_fig_AIC = figure(); hold on; 
freqs_fig_AIC.WindowState = 'maximized';
plot_histo(ALL_EST_FREQS.AIC.LL, [0.8500 0.3250 0.0980], 'AIC', ALL_EST_FREQS.KMEANS.SAMP, [0 0.4470 0.7410], 'SAMP (ours)' , BW_freqs, params.true_freqs)

% subplot(2,2,3); hold on
freqs_fig_SV_G = figure(); hold on; 
freqs_fig_SV_G.WindowState = 'maximized';
plot_histo(ALL_EST_FREQS.GAP.SV, [0.9290 0.6940 0.1250], 'GAP', ALL_EST_FREQS.KMEANS.SAMP, [0 0.4470 0.7410], 'SAMP (ours)' , BW_freqs, params.true_freqs)

% subplot(2,2,4); hold on
freqs_fig_SV_E = figure(); hold on; 
freqs_fig_SV_E.WindowState = 'maximized';
plot_histo(ALL_EST_FREQS.EFF.SV, [0.4940 0.1840 0.5560], 'EFF', ALL_EST_FREQS.KMEANS.SAMP, [0 0.4470 0.7410], 'SAMP (ours)' , BW_freqs, params.true_freqs)

hold off; 
end

end

%% Aux Function

function plot_histo(data1, color1, name1, data2, color2, name2, BW, true_params)
h1 = histogram( (data1),'BinWidth', BW,'Normalization' , 'pdf', 'FaceColor',  color1, 'EdgeColor', 'k','FaceAlpha', 1);
h2 = histogram( (data2), 'BinWidth', BW,'Normalization' , 'pdf', 'FaceColor', color2, 'EdgeColor', 'k','FaceAlpha', 1);

XI = min(true_params)-1:BW:max(true_params)+1;
num_of_comp = length(true_params);
gm1 = fitgmdist(rmoutliers(data1, 'percentiles',[2 98]), num_of_comp, 'CovarianceType', 'diagonal');   
gm2 = fitgmdist(rmoutliers(data2, 'percentiles',[2 98]), num_of_comp, 'CovarianceType', 'diagonal');

pdf1 = pdf(gm1, XI');
pdf2 = pdf(gm2, XI');

% % Plot the pdf's peaks 
% [~, pdf1_max_ind] = findpeaks(pdf1);
% [~, pdf2_max_ind] = findpeaks(pdf2);
% xline(XI(pdf1_max_ind), '--k', 'LineWidth', 4, 'color', color1) ; 
% xline(XI(pdf2_max_ind), '--k', 'LineWidth', 4, 'color', color2); 

h = plot(XI, pdf1, '-', 'Color', 'black'); % Thicker line for edge
% Overlay the main line (colored line)
LW = get(h, 'LineWidth');
h = plot(XI, pdf1, '-', 'Color', color1, 'LineWidth', LW/2); % Thinner line for the main color
set(h, 'markerfacecolor', get(h, 'color'));


h = plot(XI, pdf2, '-', 'Color', 'black'); % Thicker line for edge
% Overlay the main line (colored line)
LW = get(h, 'LineWidth');
h = plot(XI, pdf2, '-', 'Color', color2, 'LineWidth', LW/2); % Thinner line for the main color
set(h, 'markerfacecolor', get(h, 'color'));

% Plot theta_1 and theta_2 for reference
for i=1:length(true_params)
    xline(true_params(i), '--k', 'LineWidth', 4); 
end

% xlabel('$\theta$ [rad/s]', 'Interpreter','latex');
if strcmp(name1, 'GAP')
    legend([h1, h2], {name1, name2}, 'Location', 'best');
else
    legend([h1, h2], {name1, name2}, 'Location', 'north');
end

% set the figure axis
xlim([min(true_params) - 0.25, max(true_params) + 0.25])
y_lim_max = max(max(6, max(h1.Values)), max(h2.Values) );
ylim([0 y_lim_max])
yticks(0:5:y_lim_max)
xticklabels([])
xticks([min(true_params) max(true_params)])
a = get(gca,'XTick');
set(gca,'XTick',a)
xticklabels({'\theta_1','\theta_2'})
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a);

end


function  [est_signal, recovery_error] = get_signal_recovery_error(amplitudes, poles, params, dist_type)

    T  = params.interval_length;
    dt = params.dt;
    est_signal = cell(size(amplitudes));
    recovery_error = zeros(size(amplitudes));
    true_signal = params.signal-params.noise;
    % sig = zeros
    for n=1:size(amplitudes,1)
        for m=1:size(amplitudes,2)
            est_signal{n,m} = sum( (amplitudes{n,m}).*poles{n,m}.^(0:T/dt),1) ;

            if strcmp(dist_type, 'l2')
                recovery_error(n,m) = (norm(est_signal{n,m} - true_signal, 2) / norm(true_signal,2))^2 ;
            elseif strcmp(dist_type, 'l1')
                recovery_error(n,m) = norm(est_signal{n,m} - true_signal, 1) / norm(true_signal,1);
            else
                error('Unsupported distance')
            end
        end
    end

end

function [pole_estim_err, freqs_est, frequency_estim_err, damps_est, damping_estim_err]  = get_best_params_error(poles_est, params, param_num, dist_type, x_axis_param) 
    
    
    dt = params.dt;
    pole_estim_err      = zeros(size(poles_est));
    frequency_estim_err = zeros(size(poles_est));
    damping_estim_err   = zeros(size(poles_est));
    freqs_est = cell(size(poles_est));
    damps_est = cell(size(poles_est));
    
   
    % sig = zeros
    for n=1:size(poles_est,1)
        for m=1:size(poles_est,2)

            if strcmp(x_axis_param, 'FREQS_DIST')
                N = 1 +params.interval_length / params.dt;
                curr_freqs = [min(params.true_freqs) + params.FREQ_DIST_RANGE(n) * pi / (N*params.dt),min(params.true_freqs)];
                true_frequency = curr_freqs(param_num);
                true_damping   = params.true_damps(param_num);

            elseif strcmp(x_axis_param, 'SAMPLES')
                N = 1+ params.INTERVAL_LENGTH_RANGE(n)/params.dt;
                curr_freqs = [min(params.true_freqs) + 1.5 * pi / (N*params.dt), min(params.true_freqs)];
                true_frequency = curr_freqs(param_num);
                true_damping   = params.true_damps(param_num);

            elseif strcmp(x_axis_param, 'SNR')
                true_frequency = params.true_freqs(param_num);
                true_damping   = params.true_damps(param_num);
            end
        
            lambda_true = exp(-true_damping*dt + 1i*true_frequency*dt);

            freqs_est{n,m} = imag(log(poles_est{n,m}))./dt;
            damps_est{n,m} = real(log(poles_est{n,m}))./dt;

            
            % The {n,m} cell can have K frequencies / dampings . We wish to compare it
            % only to the "true_frequency", hance we take the closest
            % recovered frequencies out of those K frequencies. We choose it by the closest pole
            [~,I] = min(abs(poles_est{n,m} - lambda_true));


            if strcmp(dist_type, 'rmse')
                pole_estim_err(n,m)      = (   abs(poles_est{n,m}(I) - lambda_true)      )^2;
                frequency_estim_err(n,m) = (   abs(freqs_est{n,m}(I) - true_frequency)   )^2;
                damping_estim_err(n,m)   = (   abs(damps_est{n,m}(I) - (-true_damping))  )^2;
            elseif strcmp(dist_type, 'var')
                pole_estim_err(n,m)      = (   poles_est{n,m}(I) );
                frequency_estim_err(n,m) = (   freqs_est{n,m}(I) );
                damping_estim_err(n,m)   = (   damps_est{n,m}(I) ); 

            elseif strcmp(dist_type, 'bias')
            pole_estim_err(n,m)      = (   (poles_est{n,m}(I) - lambda_true)             );
            frequency_estim_err(n,m) = (   (freqs_est{n,m}(I) - true_frequency)          );
            damping_estim_err(n,m)   = (   (damps_est{n,m}(I) - (-true_damping))         );

            else
                error('Unsupported distance')
            end

        end
    end
end

function fig = produce_single_param_plots(Single_param_est, Single_param_bias ,crb, params, comp_num, error_type, param_name, X_AXIS)

    %set the correct name and y_limit for the figure
    if strcmp(param_name, 'freqs')
        Error_name = '$\hat{\theta}$';
    elseif strcmp(param_name, 'damps')
        Error_name = '$\hat{\alpha}$';
    else
        error('unsupported parameter name')
    end
    
    % for som data we dont need log scale
    if strcmp(error_type, 'bias')
        y_axis_conversion = @(y) y ; % no need in log scale for bias
        lgd_location = 'northeast';
    else
        y_axis_conversion = @(y) 10*log10(y);
        lgd_location = 'northeast';
    end
    
    set_IEEE_TSP_default_figure_settings()
       
    fig = figure(); hold on; 
    fig.WindowState = 'maximized';
    
    h_MDL = plot(X_AXIS.x, y_axis_conversion(Single_param_est.MDL.LL) ,'p-', 'color', [0.4660 0.6740 0.1880]);
    set(h_MDL, 'markerfacecolor', get(h_MDL, 'color'));
    
    h_AIC = plot(X_AXIS.x, y_axis_conversion(Single_param_est.AIC.LL) , 's-','color', [0.8500 0.3250 0.0980]);
    set(h_AIC, 'markerfacecolor', get(h_AIC, 'color'));
    
    h_SV_G = plot(X_AXIS.x, y_axis_conversion(Single_param_est.GAP.SV) ,'o-', 'color', [0.9290 0.6940 0.1250]);
    set(h_SV_G, 'markerfacecolor', get(h_SV_G, 'color'));
                    
    h_SV_E = plot(X_AXIS.x, y_axis_conversion(Single_param_est.EFF.SV) ,'o-', 'color', [0.4940 0.1840 0.5560]);
    set(h_SV_E, 'markerfacecolor', get(h_SV_E, 'color'));
                       
    % h = plot(X_AXIS.x, y_axis_conversion(Single_param_est.KMEANS.SV) ,'o-', 'color', [0.6350 0.0780 0.1840]);
    % set(h, 'markerfacecolor', get(h, 'color'));

    % h_fft=plot(X_AXIS.x, y_axis_conversion(Single_param_est.GAP.SAMP) ,'^-', 'color',[0.3010 0.7450 0.9330]);
    % set(h_fft, 'markerfacecolor', get(h_fft, 'color'));

    h_MF = plot(X_AXIS.x, y_axis_conversion(Single_param_est.KMEANS.SAMP) ,'-^', 'color', [0 0.4470 0.7410] );
    set(h_MF, 'markerfacecolor', get(h_MF, 'color'));

    % Add the freq diff
    plot_freqs_diff = strcmp(param_name, 'freqs') & strcmp(error_type, 'RMSE') & strcmp(X_AXIS.x_axis_param, 'SNR');
    if plot_freqs_diff
        h_freqs_diff = plot(X_AXIS.x, 10*log10( ones(size(X_AXIS.x))*abs(0.5*diff(params.true_freqs)) ) ,'color', [0.5 0.5 0.5], 'LineStyle', '--');
    end
    
    % add CRB for variance or RMSE
    if ~isempty(crb)
        if isempty(Single_param_bias)
            bias = 0;
        else
            bias = Single_param_bias.GAP.SAMP;
        end
        if strcmp(error_type, 'RMSE')
            h_CRB=plot(X_AXIS.x, 10*log10(sqrt(crb + bias.^2) ) ,'color', 'black');
        elseif strcmp(error_type, 'var')
            h_CRB=plot(X_AXIS.x, 10*log10(crb) ,'color', 'black');
        end
        set(h_CRB, 'markerfacecolor', get(h_CRB, 'color'));
    end

    % Set the legend
    if ~isempty(crb) & plot_freqs_diff
        lgd = legend([h_CRB, h_freqs_diff, h_MDL, h_AIC, h_SV_G, h_SV_E, h_MF], {'CRB','$|\theta_1-\theta_2|/2$', 'MDL', 'AIC', 'GAP', 'EFF', 'SAMP (ours)'}, 'Location', lgd_location, 'Interpreter','latex');    
    elseif  isempty(crb) & plot_freqs_diff
        lgd = legend([h_freqs_diff, h_MDL, h_AIC, h_SV_G, h_SV_E, h_MF], {'$|\theta_1-\theta_2|/2$', 'MDL', 'AIC', 'GAP', 'EFF', 'SAMP (ours)'}, 'Location', lgd_location, 'Interpreter','latex');    
    elseif  ~isempty(crb) & ~plot_freqs_diff
        lgd = legend([h_CRB, h_MDL, h_AIC, h_SV_G, h_SV_E, h_MF], {'CRB', 'MDL', 'AIC', 'GAP', 'EFF', 'SAMP (ours)'}, 'Location', lgd_location, 'Interpreter','latex');    
   
    else
        lgd = legend([h_MDL, h_AIC, h_SV_G, h_SV_E, h_MF], {'MDL', 'AIC', 'GAP', 'EFF', 'SAMP (ours)'}, 'Location', lgd_location, 'Interpreter','latex');   
    end
    
    % set log scale
    % if strcmp(X_AXIS.x_axis_param, 'SNR') & strcmp(error_type, 'RMSE')
    %     set(gca,'YScale','log');
    % end

    % Adjust fonts and axis labels
    fontsize(lgd,32,'points')      
    xlim(X_AXIS.x_axis_range)  
    xlabel(X_AXIS.x_axis_label,  'Interpreter','latex')
   
    if strcmp(error_type, 'bias')
        ylabel([error_type, '(', Error_name,')'], 'Interpreter','latex')
    else
        ylabel([error_type, '(', Error_name,')', ' [dB]'], 'Interpreter','latex')
    end
    grid on;
    box on;

    if strcmp(X_AXIS.x_axis_param, 'FREQS_DIST')
        xticks([1,2,3,4])
        xticklabels({'$\frac{\pi}{N}$', '$\frac{2\pi}{N}$', '$\frac{3\pi}{N}$', '$\frac{4\pi}{N}$'})
        set(gca, 'TickLabelInterpreter', 'latex')
    end
    % Set figure limits
    % if strcmp(error_type, 'bias')
    %     ylim([0 1]);
    %     yticks([0, 0.5, 1])
    % elseif strcmp(error_type, 'RMSE')
    %     ylim([-25 -5]);
    %     yticks([-25 -20 -15, -10, -5 ]);
    % end

end




