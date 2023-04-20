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
% Fill in data_model fields used in all analyses
data_model.nrois    = nrois;    % Number of ROIs
data_model.nruns    = nruns;    % Number of iterations of SIR dynamics 
data_model.gammas   = gammas;   % Array of values of the recovery probability
data_model.rhos     = rhos;     % Array of values of the network density
data_model.nreps    = nreps;    % Number of repetitions
data_model.pout_dyn = pout_dyn; % Output path for SIR data 