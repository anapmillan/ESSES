function d2soz = VR_d2soz(nodes, aij, esoz_map)
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
% function d2soz = VR_d2soz(nodes, aij, esoz_map)
% Perform virtual resection of the nodes in "nodes" and return the seed
% efficiency.
%   INPUTS:
%       nodes:      Set of node indices to be removed.
%       aij:        Adjacency matrix (unresected).
%       esoz_map:   Seed probability map (to derive the seed).
%   OUTPUTS:
%       d2soz:      Seed efficiency.

    aij_VR = prune_nodes(aij, nodes);
    d2soz = matrix2distance(aij_VR, esoz_map);