function [minimum_solution, optimal_solution] = ...
    simulated_annealing(handle_neighbourhood, handle_loss,...
    initial_connection_indices, data_SA)
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
%function [minimum_solution, optimal_solution] = ...
%    simulated_annealing(handle_neighbourhood, handle_loss,...
%    initial_connection_indices, data) 
% Script to run simulated annealing
% The variable to be optimised is the list of connections to be removed: 
% connection_indices.
% The optimization function is given by the cost function.


    options = data_SA.options;
    options.Generator = handle_neighbourhood;


%    options.number_of_connections_in_list = size(out_soz_edges,1);
%     options.InitTemp = 0.1;
%     options.StopTemp = 1e-3;
%     options.StopVal = 0;
%     options.Verbosity = 2;
%     options.CoolSched = @(T) (.6*T);
%     options.MaxConsRej = 500;

    % Run the SA
    loss = handle_loss;
    [minimum_solution, optimal_solution] = ...
        anneal(loss, initial_connection_indices, options);