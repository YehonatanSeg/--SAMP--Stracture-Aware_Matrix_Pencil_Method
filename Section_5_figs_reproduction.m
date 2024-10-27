%% General Settings
Nreps = 500;
damps = [0,0] ;
noise_distribution = 'normal'; % 'bi_normal', 't_student', 'uniform';
addpath('toolbox')
%% Generate Figs 3.a, 3.b

% Generate Data
% For 3.a use damps = [0,0], for 3.b use damps = [0.03,0.05];
[Poles, Left_modes,  MP_Amplitudes, NEW_Amplitudes,...
    SDD, GAP, EFF, KMEANS, MDL, AIC, CRB, params] = gen_data_SNR(noise_distribution, Nreps, damps);

% Preprocess data
[SDD, GAP, KMEANS] = preprocess_data(params, Poles, Left_modes, SDD, GAP, KMEANS);

% Produce graph
plot_detection_graph(params, AIC, MDL, GAP,EFF, KMEANS, 'SNR')

% Force the plot to update
drawnow; 
pause(0.1);
%% Generate Figs 3.c, 3.d

% Generate Data
% For 3.c use damps = [0,0], for 3.d use damps = [0.03,0.05];
[Poles, Left_modes,  MP_Amplitudes, NEW_Amplitudes,...
    SDD, GAP, EFF, KMEANS, MDL, AIC, CRB, params] = gen_data_N(noise_distribution, Nreps, damps);

% Preprocess data
[SDD, GAP, KMEANS] = preprocess_data(params, Poles, Left_modes, SDD, GAP, KMEANS);

% Produce graph
plot_detection_graph(params, AIC, MDL, GAP,EFF, KMEANS, 'SAMPLES')

% Force the plot to update
drawnow; 
pause(0.1);
%% Generate Figs 3.e, 3.f

% Generate Data
% For 3.e use damps = [0,0], for 3.f use damps = [0.03,0.05];
[Poles, Left_modes,  MP_Amplitudes, NEW_Amplitudes,...
    SDD, GAP, EFF, KMEANS, MDL, AIC, CRB, params] = gen_data_FREQS_DIST(noise_distribution, Nreps, damps);

% Preprocess data
[SDD, GAP, KMEANS] = preprocess_data(params, Poles, Left_modes, SDD, GAP, KMEANS);

% Produce graph
plot_detection_graph(params, AIC, MDL, GAP,EFF, KMEANS, 'FREQS_DIST')

% Force the plot to update
drawnow; 
pause(0.1);
%% Generate Figs 4,5,8,9

% Generate Data 
% For the undamped case use damps = [0,0], forthe damped case use damps = [0.03,0.05];
[Poles, Left_modes,  MP_Amplitudes, NEW_Amplitudes,...
    SDD, GAP, EFF, KMEANS, MDL, AIC, CRB, params] = gen_data_SNR(noise_distribution, Nreps, damps);

% Preprocess data
[SDD, GAP, KMEANS] = preprocess_data(params, Poles, Left_modes, SDD, GAP, KMEANS);

% Produce graph
plot_estimation_graphs(params, CRB, AIC, MDL, GAP,EFF, KMEANS, 'SNR')

%% Generate Figs 10.a, 10.b, 10.c, 10.d

% Generate Data 
% For 10.a and 10.b use damps = [0,0], for 10.c and 10.d use damps = [0.03,0.05];
[Poles, Left_modes,  MP_Amplitudes, NEW_Amplitudes,...
    SDD, GAP, EFF, KMEANS, MDL, AIC, CRB, params] = gen_data_N(noise_distribution, Nreps, damps);

% Preprocess data
[SDD, GAP, KMEANS] = preprocess_data(params, Poles, Left_modes, SDD, GAP, KMEANS);

% Produce graph
plot_estimation_graphs(params, CRB, AIC, MDL, GAP,EFF, KMEANS, 'SAMPLES')

% Force the plot to update
drawnow; 
pause(0.1);
%% Generate Figs 11.a, 11.b, Figs 11.c, 11.d

% Generate Data 
% For 11.a and 11.b use damps = [0,0], for 11.c and 11.d use damps = [0.03,0.05];
[Poles, Left_modes,  MP_Amplitudes, NEW_Amplitudes,...
    SDD, GAP, EFF, KMEANS, MDL, AIC, CRB, params] = gen_data_FREQS_DIST(noise_distribution, Nreps, damps);

% Preprocess data
[SDD, GAP, KMEANS] = preprocess_data(params, Poles, Left_modes, SDD, GAP, KMEANS);

% Produce graph
plot_estimation_graphs(params, CRB, AIC, MDL, GAP,EFF, KMEANS, 'FREQS_DIST')

% Force the plot to update
drawnow; 
pause(0.1);
%% Generate Fig 6.a, 6.b

% Generate Data 
% For 6.a use damps = [0,0], for 6.b use damps = [0.03,0.05];
[Poles, Times,NOF, params] = gen_data_time_comparisons_MDL_AIC_SAMP(noise_distribution, Nreps, damps);

% Produce graph
plot_time_compare_graphs_MDL_AIC_SAMP(Poles, Times, params)

% Force the plot to update
drawnow; 
pause(0.1);
%% Generate Fig 7a, 7.b

% Generate Data 
% For 7.a use damps = [0,0], for 7.b use damps = [0.03,0.05];
[params, TIMES, EST] = gen_data_amps_comparisons_MP_SAMP(noise_distribution, Nreps, damps);

% Produce graph
plot_amps_compare_graphs_MP_SAMP(params, TIMES, EST)