function [wij, seed, pattern, castor_map] = ...
    prepare_data_patient(pat_dir, get_SEEG, get_map)
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
% function [wij, seed, pattern] = prepare_data_patient(pat_dir,get_SEEG, get_map)
%   Read patient data. 
%   INPUTS:
%       pat_dir: struct variable indicating the filenames of the patient
%                files. Field options include: 
%                seed_filename:    txt file including the seed regions as a
%                                  column vector.
%                network_filename: txt file including the network matrix as 
%                                  an [Nr_ROIs x Nr_ROIs] matrix 
%                pattern_filename: mat file with the iEEG spreading pattern
%                map_filename:     txt file including the seed_probability
%                                  map as an [Nr_ROIs x 1] vector.
%       (optional) get_SEEG, get_map: tags to control whether the iEEG and
%                seed_probability data are retrieved (default = 0).
%  OUTPUTS:
%      wij:         Patient-specific brain connectivity, [Nr_ROIs x
%                   Nr_ROIs] matrix.
%      seed:        Resection area.
%      pattern:     Structure variable enconding the iEEG seizure propagation
%                   pattern.
%      castor_map:  [Nr_ROIs x 1] vector indicating the probability that
%                   each ROI initiates a seizure.

    
    if nargin<2
        get_map = 0;
        if nargin<1
            get_SEEG=0;
        end
    end

    % Read seed file:
    seed    = dlmread(pat_dir.seed_filename);
    % Read network file:
    wij     = dlmread(pat_dir.network_filename);
    wij     = (wij+wij')/2; %make sure that the network is symmetric

    % iEEG seizure propagation pattern:
    if get_SEEG
        load(pat_dir.pattern_filename,'pattern');
        pattern = pattern(end);
    else
        pattern = [];
    end

    % Seed-probability map:
    if get_map
        castor_map = dlmread(pat_dir.map_filename);
    else
        castor_map = [];
    end