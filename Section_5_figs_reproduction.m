%% General Settings
Nreps = 500;
noise_distribution = 'normal'; % 'bi_normal', 't_student', 'uniform';
addpath('toolbox')
%% Generate Figs 3.a, 3.b

% For 3.a use damps = [0,0], for 3.b use damps = [0.03,0.05];
damps = [0,0] ;

% Generate Data 
[Poles, Left_modes,  MP_Amplitudes, NEW_Amplitudes,...
    SDD, GAP, EFF, KMEANS, MDL, AIC, CRB, params] = gen_data_SNR(noise_distribution, Nreps, damps);

% Preprocess data
[SDD, GAP, KMEANS] = preprocess_data(params, Poles, Left_modes, SDD, GAP, KMEANS);

% Produce graph
plot_detection_graph(params, AIC, MDL, GAP,EFF, KMEANS, 'SNR')

%% Generate Figs 3.c, 3.d

% For 3.c use damps = [0,0], for 3.d use damps = [0.03,0.05];
damps = [0,0] ;

% Generate Data 
[Poles, Left_modes,  MP_Amplitudes, NEW_Amplitudes,...
    SDD, GAP, EFF, KMEANS, MDL, AIC, CRB, params] = gen_data_N(noise_distribution, Nreps, damps);

% Preprocess data
[SDD, GAP, KMEANS] = preprocess_data(params, Poles, Left_modes, SDD, GAP, KMEANS);

% Produce graph
plot_detection_graph(params, AIC, MDL, GAP,EFF, KMEANS, 'SAMPLES')

%% Generate Figs 3.e, 3.f

% For 3.e use damps = [0,0], for 3.f use damps = [0.03,0.05];
damps = [0,0] ;

% Generate Data 
[Poles, Left_modes,  MP_Amplitudes, NEW_Amplitudes,...
    SDD, GAP, EFF, KMEANS, MDL, AIC, CRB, params] = gen_data_FREQS_DIST(noise_distribution, Nreps, damps);

% Preprocess data
[SDD, GAP, KMEANS] = preprocess_data(params, Poles, Left_modes, SDD, GAP, KMEANS);

% Produce graph
plot_detection_graph(params, AIC, MDL, GAP,EFF, KMEANS, 'FREQS_DIST')


%% Generate Figs 4,5,8,9

% For the undamped case use damps = [0,0], forthe damped case use damps = [0.03,0.05];
damps = [0,0] ;

% Generate Data 
damps = [0,0]
[Poles, Left_modes,  MP_Amplitudes, NEW_Amplitudes,...
    SDD, GAP, EFF, KMEANS, MDL, AIC, CRB, params] = gen_data_SNR(noise_distribution, Nreps, damps);

% Preprocess data
[SDD, GAP, KMEANS] = preprocess_data(params, Poles, Left_modes, SDD, GAP, KMEANS);

% Produce graph
plot_estimation_graphs(params, CRB, AIC, MDL, GAP,EFF, KMEANS, 'SNR')

%% Generate Figs 10.a, 10.b, 10.c, 10.d

% For 10.a and 10.b use damps = [0,0], for 10.c and 10.d use damps = [0.03,0.05];
damps = [0,0] ;

% Generate Data 
[Poles, Left_modes,  MP_Amplitudes, NEW_Amplitudes,...
    SDD, GAP, EFF, KMEANS, MDL, AIC, CRB, params] = gen_data_N(noise_distribution, Nreps, damps);

% Preprocess data
[SDD, GAP, KMEANS] = preprocess_data(params, Poles, Left_modes, SDD, GAP, KMEANS);

% Produce graph
plot_estimation_graphs(params, CRB, AIC, MDL, GAP,EFF, KMEANS, 'SAMPLES')


%% Generate Figs 11.a, 11.b, Figs 11.c, 11.d

% For 11.a and 11.b use damps = [0,0], for 11.c and 11.d use damps = [0.03,0.05];
damps = [0,0] ;

% Generate Data 
[Poles, Left_modes,  MP_Amplitudes, NEW_Amplitudes,...
    SDD, GAP, EFF, KMEANS, MDL, AIC, CRB, params] = gen_data_FREQS_DIST(noise_distribution, Nreps, damps);

% Preprocess data
[SDD, GAP, KMEANS] = preprocess_data(params, Poles, Left_modes, SDD, GAP, KMEANS);

% Produce graph
plot_estimation_graphs(params, CRB, AIC, MDL, GAP,EFF, KMEANS, 'FREQS_DIST')

%% Generate Fig 6.a, 6.b

% For 6.a use damps = [0,0], for 6.b use damps = [0.03,0.05];
damps = [0,0] ;

% Generate Data 
[Poles, Times,NOF, params] = gen_data_time_comparisons_MDL_AIC_SAMP(noise_distribution, Nreps, damps);

% Produce graph
plot_time_compare_graphs_MDL_AIC_SAMP(Poles, Times, params)

%% Generate Fig 7a, 7.b

% For 7.a use damps = [0,0], for 7.b use damps = [0.03,0.05];
damps = [0,0] ;

% Generate Data 
 [params, TIMES, EST] = gen_data_amps_comparisons_MP_SAMP(noise_distribution, Nreps, damps);

% Produce graph
plot_amps_compare_graphs_MP_SAMP(params, TIMES, EST)