function [SDD, GAP, KMEANS] = preprocess_data(params, Poles, Left_modes, SDD, GAP, KMEANS)

    m=1;
    for i=1:size(Left_modes,1)
        disp([num2str(round(100*m / (size(Left_modes,1)*size(Left_modes,2)))) '%'])
        for j=1:size(Left_modes,2)
           
                field = 'SAMP';

                % Extract the relevant information from the left mode
                epsilon = compute_signal_detection_features(Left_modes{i,j}, Poles{i,j}, params);
                epsilon    = epsilon/max(epsilon);
                
                % Get the indices according to each method
                ind_SDD                         = get_true_idx(epsilon, 'SD');
                ind_GAP                         = get_true_idx(epsilon, 'GAP');
                ind_KMEANS                      = get_true_idx(epsilon, 'kmeans');

                % Poles estimation
                SDD.Poles.(field){i,j}          = Poles{i,j}(ind_SDD);
                GAP.Poles.(field){i,j}          = Poles{i,j}(ind_GAP);
                KMEANS.Poles.(field){i,j}       = Poles{i,j}(ind_KMEANS);

                % Number of components (detection)
                SDD.NOF.(field)(i,j)            = length(ind_SDD);
                GAP.NOF.(field)(i,j)            = length(ind_GAP);
                KMEANS.NOF.(field)(i,j)         = length(ind_KMEANS);

                m=m+1;
        end
    end
end


