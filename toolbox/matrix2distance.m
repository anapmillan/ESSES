function rmin = matrix2distance(aij, esoz_map)
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
% function rmin = matrix2distance(aij, esoz_map)
%   Measure seed efficiency. 
%   INPUTS:
%       aij:        Adjency matrix (target, either original or resected).
%       esoz_map:   Seed-probability map.
%   OUTPUT:
%       rmin:       Seed-efficiency.
% The seed is defined as the non-zero indices of esoz_map.
 
    % Seed regions:
    esoz0           = find(esoz_map>0);
    
    % Non-seed nodes:
    ne_soz          = 1:size(aij,1);
    ne_soz(esoz0)  = [];
    % Distances should not include seed nodes for stability

    % Graph object:
    gij     = graph(aij);
    % Distance matrix
    dij      = 1./distances(gij, esoz0);
    % Keep only non-seed nodes
    dij      = dij(:,ne_soz);
    esoz_w   = esoz_map(esoz0); %lower efficiency for lower charge
    dijw     = dij.*esoz_w;
    rmin     = mean(dijw(:));

 


