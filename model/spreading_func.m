function [order_in, order_out, frac_X] =...
    spreading_func(aij, data_run, seed)
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
% function [order_in, order_out, frac_x] = ...
%          spreading_func(aij, data_run, new_seed)
% Prepare parameters and filenames to call the SIR spreading function
% ("function_SIR_adaptive_inverted.m").
% INPUTS:
%   aij:       [Nr_ROIs x Nr_ROIs] weighted and symmetric adjacency matrix
%   data_run:  Data for SIR simulations. Required files:
%       data_run.is_VR:  tag for VR analyses, used in the filenames 
%                        (e.g. 'BS','VR','') 
%       data_run.km0:    Network mean degree (for filenames) 
%       data_run.gamma:  Recovery probability
%       data_run.rep:    Index of realization (for filenames)
%   new_seed:            Spreading seed.
% IMPLICIT INPUTS (global variables):
%   data_model.nruns:    Number of SIR realizations
%   data_model.pout_dyn: Path to save SIR data.
% OUTPUS:
%   order_in:  [nr_ROIs x nr_ROIs+1] matrix where each element
%              (i,j) indicates the probability that ROI i becomes infected
%              at step j
%   order_out: [nr_ROIs x nr_ROIs+1] matrix as order_in, but indicating the 
%              probability of recovery
%   frac_x:    [3 x nruns] matrix indicating the fraction of ROIs in the S,
%              I, and R state in the final step of the simulation.


    global data_model 
    s_file = 1;
    %%

    name_files = sprintf('%s_km%.4f_gamma%.4f_nruns%d_rep%d',...
        data_run.is_VR, data_run.km0, data_run.gamma,...
        data_model.nruns, data_run.rep);

    n_out   = sprintf('%s/dyn%s_%s.txt', ...
        data_model.pout_dyn, data_run.is_VR,  name_files);

    [order_in, order_out, frac_X] = func_SIR_adaptive_inverted(...
        aij, data_run.gamma, seed);
 
    order_in = order_in';
    order_out = order_out';

    if s_file; dlmwrite(n_out, order_in, ' '); end