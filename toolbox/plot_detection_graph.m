function plot_detection_graph(params, AIC, MDL, GAP,EFF, KMEANS, x_axis_param)

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
    X_AXIS.x_axis_label = '$|\theta_1-\theta_2|$ [rad/s]';
elseif strcmp(X_AXIS.x_axis_param, 'SAMPLES')
    X_AXIS.x = params.INTERVAL_LENGTH_RANGE / params.dt;
    X_AXIS.x_axis_range = [X_AXIS.x(1), X_AXIS.x(end)];
    X_AXIS.x_axis_label = 'N';
end

%% Create the detection probability of each method
if strcmp(X_AXIS.x_axis_param, 'NUM') 
    M = params.NUM_FREQS_RANGE';
else
    M = length(params.true_freqs); 
end
Pd = [];
s_ind = find(X_AXIS.x == X_AXIS.x_axis_range(1));
e_ind = find(X_AXIS.x == X_AXIS.x_axis_range(end));
for k= 1:numel(METHODS)
    method = eval(METHODS{k});
    fields = fieldnames(method.Poles);
    for i =1:length(fields)
        Pd.(METHODS{k}).(fields{i})  = mean(method.NOF.(fields{i}) == M, 2);
        AUC.(METHODS{k}).(fields{i}) = sum(Pd.(METHODS{k}).(fields{i})(s_ind: e_ind)) / length(X_AXIS.x(s_ind: e_ind));
    end                                            
end


%%
fig = figure(); hold on;
fig.WindowState = 'maximized';

h = plot(X_AXIS.x, Pd.MDL.LL ,'p-', 'color', [0.4660 0.6740 0.1880]);
set(h, 'markerfacecolor', get(h, 'color')); 

h = plot(X_AXIS.x, Pd.AIC.LL , 's-','color', [0.8500 0.3250 0.0980]);
set(h, 'markerfacecolor', get(h, 'color')); 

h = plot(X_AXIS.x, Pd.GAP.SV ,'o-', 'color', [0.9290 0.6940 0.1250]);
set(h, 'markerfacecolor', get(h, 'color'));

h = plot(X_AXIS.x, Pd.EFF.SV ,'o-', 'color', [0.4940 0.1840 0.5560]);
set(h, 'markerfacecolor', get(h, 'color'));

h = plot(X_AXIS.x, Pd.KMEANS.SAMP ,'-^', 'color', [0 0.4470 0.7410]);
set(h, 'markerfacecolor', get(h, 'color'));

lgd = legend(['MDL                ',num2str(sprintf('%.2f',AUC.MDL.LL))],...
             ['AIC                  ',num2str(sprintf('%.2f',AUC.AIC.LL))],...
             ['GAP                 ',num2str(sprintf('%.2f',AUC.GAP.SV))],...
             ['EFF                  ',num2str(sprintf('%.2f',AUC.EFF.SV))],...
             ['SAMP (our)     ',num2str(sprintf('%.2f',AUC.KMEANS.SAMP))],...
             'Location','northeast');

fontsize(lgd,40,'points')      
xlim(X_AXIS.x_axis_range)  
xlabel(X_AXIS.x_axis_label,  'Interpreter','latex')
ylabel('$p_d$',  'Interpreter','latex')

if strcmp(X_AXIS.x_axis_param, 'DIFF')
    xticks([1,2,3,4])
    xticklabels({'$\frac{\pi}{N}$', '$\frac{2\pi}{N}$', '$\frac{3\pi}{N}$', '$\frac{4\pi}{N}$'})
    set(gca, 'TickLabelInterpreter', 'latex')
end

end

