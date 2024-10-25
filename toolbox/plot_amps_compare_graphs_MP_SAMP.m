function plot_amps_compare_graphs_MP_SAMP(params, TIMES, EST)
% (params, new_amps_times,MP_amps_times, new_amps_est1,new_amps_est2, MP_amps_est1,MP_amps_est2)
    
    set_IEEE_TSP_default_figure_settings()
    fig = figure(); hold on
    fig.WindowState = 'maximized';
    set(gcf, 'Units','centimeters')    

    %% Averaging resutls
    % Times
    SAMP_amps_times_mean   = mean(TIMES.SAMP ,2);
    MP_amps_times_mean     = mean(TIMES.MP,2);

    % Amps
    SAMP_amps_RMSE         = rmse(EST.SAMP_est1, params.true_amps(1),2) + rmse(EST.SAMP_est1, params.true_amps(2),2);
    MP_amps_RMSE           = rmse(EST.MP_est1,  params.true_amps(1),2) + rmse(EST.MP_est2,  params.true_amps(2),2);

    h=plot(MP_amps_times_mean, 10*log10(MP_amps_RMSE), 's-', 'color', [0.8500 0.3250 0.0980]); 
    set(h, 'markerfacecolor', get(h, 'color'));

    h=plot(SAMP_amps_times_mean, 10*log10(SAMP_amps_RMSE),'^-', 'color', [0 0.4470 0.7410]); 
    set(h, 'markerfacecolor', get(h, 'color'));

    lgd = legend('MP', 'SAMP (ours)');
    % fontsize(lgd,40,'points')    
    xlabel('Computation time [s]')
    ylabel('RMSE [dB]')
    yticks([-12, -7, -2])
    ylim([-12, -1.5])
    % set(gca,'YScale','log');
    % axis('tight')
    grid on;
end