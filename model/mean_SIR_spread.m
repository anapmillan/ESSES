function tpt_acu_m = mean_SIR_spread(aij, data_run, new_seed)
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
% function tpt_acu_m = mean_SIR_spread(aij, data_run, new_seed)
%   Run SIR dynamics and meausure outputs for VR.
%   INPUTS:
%       aij:        [nr_ROIs x nr_ROIs] adjacency matrix. 
%                   Weighted and symmentric.
%       data_run:   struct variable with data for spreading_func.m
%       new_seed:   seed of the SIR dynamics.
%  OUTPUS: 
%       tpt_acu_m:  Fraction of nodes that became infected (I+R).

    [~, ~, frac_x] = spreading_func(aij, data_run, new_seed);
    % fracion of nodes in each state, on each run, final step
    % susteptible infected recovered
    frac_x      = mean(frac_x,2);
    tpt_acu_m   = 1-frac_x(1);

