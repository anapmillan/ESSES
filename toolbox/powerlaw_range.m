function n_ch = powerlaw_range(xx)
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
% function n_ch = powerlaw_range(xx) 
% Generate integers according to a power-law distribution with exponent -1.
    expp    = -1;
    pi      = xx.^expp;
    pi      = pi/sum(pi);
    pi_acu  = zeros(size(pi));
    for i   = 1:numel(xx)
        pi_acu(i) = sum(pi(1:i));
    end
    aux     = rand;
    n_ch    = min(find(pi_acu>aux));