addpath(fullfile(fileparts(pwd), 'toolbox'));
clear;
close all;
%% General Settings
Nreps = 1000;
damps = [0,0,0,0] ;
noise_distribution = 'normal'; % 'bi_normal', 't_student', 'uniform';
%% Generate Figs 2.a, 2.b

% Generate Data
% For 2.a use damps = [0,0,0,0], for 2.b use damps = [0.03,0.05,0.03,0.05];
[Poles, Left_modes, Amplitudes,...
    SDD, GAP, EFF, KMEANS, MDL, AIC, MAP, EVT, params] = gen_data_SNR_SAMP_pp(noise_distribution, Nreps, damps);

% Preprocess data
[SDD, GAP, KMEANS] = preprocess_data(params, Poles, Left_modes, SDD, GAP, KMEANS);

% Produce graph
plot_detection_graph_SAMP_pp(params, EVT, MDL, GAP, SDD, KMEANS, 'SNR')

% Force the plot to update
drawnow; 
pause(1);
%% Generate Figs 2.c, 2.d

% Generate Data
% For 2.c use damps = [0,0,0,0], for 2.d use damps = [0.03,0.05,0.03,0.05];
[Poles, Left_modes, Amplitudes,...
    SDD, GAP, EFF, KMEANS, MDL, AIC, MAP, EVT, params] = gen_data_CLUS_DIST_SAMP_pp(noise_distribution, Nreps, damps);

% Preprocess data
[SDD, GAP, KMEANS] = preprocess_data(params, Poles, Left_modes, SDD, GAP, KMEANS);

% Produce graph
plot_detection_graph_SAMP_pp(params, EVT, MDL, GAP, SDD, KMEANS, 'CLUST')

% Force the plot to update
drawnow; 
pause(1);
%% Generate Figs 2.e, 2.f

% Generate Data
% For 2.e use damps = [0,0,0,0], for 2.f use damps = [0.03,0.05,0.03,0.05];
[Poles, Left_modes, Amplitudes,...
    SDD, GAP, EFF, KMEANS, MDL, AIC, MAP, EVT, params] = gen_data_N_SAMP_pp(noise_distribution, Nreps, damps);

% Preprocess data
[SDD, GAP, KMEANS] = preprocess_data(params, Poles, Left_modes, SDD, GAP, KMEANS);

% Produce graph
plot_detection_graph_SAMP_pp(params, EVT, MDL, GAP, SDD, KMEANS, 'SAMPLES')

% Force the plot to update
drawnow; 
pause(1);
