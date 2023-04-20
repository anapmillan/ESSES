function aij = threshold_matrix(wij, km0)
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
% function wij = threshold_matrix(wij, kmo).
% Threshold adjacency matrix to target mean degree.
% INPUTS: 
%       wij: original adjacency matrix.
%       kmo: target meen degree.
% OUTPUS:
%       aij: Thresholded adjacency matrix.

    wij = (wij+wij')/2;
    thw = sort(unique(wij), 'descend');
    thw = [1;thw];
    km = 0;
    i = 1;
    while km<km0
        aij = wij>thw(i);
        km = mean(sum(aij));
        i = i+1;
    end
    ith = i-1;
    aij = wij > thw(ith);
    aij = aij.*wij;