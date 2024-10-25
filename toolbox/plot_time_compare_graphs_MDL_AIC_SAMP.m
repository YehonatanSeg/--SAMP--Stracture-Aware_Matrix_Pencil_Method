function plot_time_compare_graphs_MDL_AIC_SAMP(Poles, Times, params)

fig = figure(); hold on
fig.WindowState = 'maximized';
set_IEEE_TSP_default_figure_settings()

% Times
MDL_times_mean        = mean(Times.MDL ,2);
AIC_times_mean        = mean(Times.AIC ,2);

SAMP_gap_times_mean     = mean(Times.SAMP_gap ,2);
SAMP_kmeans_times_mean  = mean(Times.SAMP_kmeans ,2);

fields = fieldnames(Poles);
for i = 1:length(fields)
    for k=1:size(Poles.MDL,1)
            
        % inf values:
        est1 = cellfun(@(x) x(1), Poles.(fields{i})(k,:));
        est2 = cellfun(@(x) x(2), Poles.(fields{i})(k,:));
        inf_ind = find(isinf(est2));
        est2(inf_ind) = est1(inf_ind);
        RMSE.(fields{i})(k)       = rmse( est1,  params.true_poles(1), 2) + rmse( est2,  params.true_poles(2),2);
    end
end

h=plot(MDL_times_mean, 10*log10(RMSE.MDL), 'p-', 'color', [0.4660 0.6740 0.1880]); 
set(h, 'markerfacecolor', get(h, 'color'));

h=plot(AIC_times_mean, 10*log10(RMSE.AIC), 's-','color', [0.8500 0.3250 0.0980]); 
set(h, 'markerfacecolor', get(h, 'color'));

h=plot(SAMP_kmeans_times_mean, 10*log10(RMSE.SAMP_kmeans),'-^', 'color', [0 0.4470 0.7410] ); 
set(h, 'markerfacecolor', get(h, 'color'));


lgd = legend('MDL', 'AIC', 'SAMP (ours)');
xlabel('Computation time [s]')
ylabel('RMSE [dB]')
axis('tight')
grid on;
