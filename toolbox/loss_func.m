function cost = loss_func(idx_nodes, aij, esoz_map)
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
%   function cost = loss_func_d2soz(idx_nodes, w_metric_cost)
%   Define cost function of simulated annealing.
%   The cost is defined as the efficiency of the seed following the
%   resection.
%   INPUTS:
%       idx_nodes:  Indices of nodes to be resected
%       aij:        Adjacency matrix
%       esoz_map:   Seed-probability map
%   OUTPUT:
%       cost:       Efficiency of the seed following resection of idx_nodes
    idx_nodes =  idx_nodes';
    cost = VR_d2soz(idx_nodes, aij, esoz_map);
