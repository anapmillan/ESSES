function [tpt_acu_m, seed_size] = run_VR(nseeds, data_run)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you use this code, then please cite:
% 1.- Ana P. Millan, et al. "Individualized epidemic spreading models 
%     predict epilepsy surgery outcomes: a pseudo-prospective study." 
%     medRxiv (2023): 2023-03. doi: https://doi.org/10.1101/2023.03.16.23287370
% 2.- Ana P. Millan, et al. "The role of epidemic spreading in seizure 
%     dynamics and epilepsy surgery." Network Neuroscience (2022): 1-55.
%     doi: https://doi.org/10.1162/netn_a_00305
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function tpt_acu_m = run_VR(ra_idx, nseeds, data_run)
% Run VR analysis for a set of nseeds different iterations of the seed
% regions. The resection area is given by data_patient.RA
% INPUTS:
%       nseeds:     Scalar, number of seed realizations
%       data_run:   Structure with the fields:
%           irep:   Scalar, index of iteration
% Implicit inputs in the global variables data_model and data_patient
% include:
%       data_patient.aij:   Original (thresholded) network, [Nr_ROIs x
%                           NR_ROIs] matrix
%       data_patient.RA:    1D array indicating the indices of the nodes in
%                           the resection
%       data_patient.esoz:  1D array indicating the seed.
%                           Can be either the seed nodes (size [1 x number_of_seeds])
%                           or the seed probability map (size [1 x
%                           Nr_ROIs]) according to data_model.w_seed
%       data_model.w_seed:  String or char array. If equal to 'map', then
%                           data_patient.esoz is used as the seed probability map
% OUTPUTS:
%       tpt_acu_m:  Matrix of size [nseeds x 2]. For each seed, spreading in
%                   the baseline (first column) and resected (second column)
%                   networks
%       seed_size:  Matrix of size [nseeds x 1] indicating the size of the
%                   seed in each iteration

    global data_model data_patient
    aij = data_patient.aij;

    % Apply resection to network
    aij_VR = prune_nodes(aij, data_patient.RA);

    % Loop over seeds
    tpt_acu_m = zeros(nseeds,2);
    seed_size = zeros(nseeds,1);
    for i = 1:nseeds
        data_run.rep = i;
        % Select new seed
        if ~strcmp(data_model.w_seed, 'map')
            new_seed = data_patient.esoz;
            s = numel(new_seed);
        else
            s = 0;
            while s==0
                new_seed = find(data_patient.esoz>rand(246,1))';
                s = numel(new_seed);
            end
        end
        seed_size(i) = s;

        % Run dynamics and measure spreading BASELINE and VR
        data_run.is_VR = 'BS'; %tag for the namefiles so that BS and VR are diff
        aux1 = mean_SIR_spread(aij, data_run, new_seed);
        data_run.is_VR = 'VR';
        aux2 = mean_SIR_spread(aij_VR, data_run, new_seed);

        tpt_acu_m(i,:) = [aux1, aux2];

    end