function bij = prune_nodes(aij, nodes)
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
% function bij = prune_nodes(aij, nodes)
%   Auxiliary function to perform a virtual resection. All connections of
%   the resected nodes are set to 0.
%   INPUTS:
%       aij:    [Nr_ROIs x Nr_ROIs] adjacency matrix
%       nodes:  indices of nodes to be resected
%   OUTPUS:
%       bij:    [Nr_ROIs x Nr_ROIs] resected adjacency matrix
    bij = aij;
    bij(nodes,:) = 0;
    bij(:,nodes) = 0;