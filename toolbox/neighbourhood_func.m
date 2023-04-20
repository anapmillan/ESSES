function new_idx = neighbourhood_func(connection_indices, data_model)
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
% function new_idx = neighbourhood_func(connection_indices, data_model) 
% Select new resection nodes
% INPUTS:
%   connection_indices:         indices of current resection
%   data_model:                 structure function including the field:
%       data_model.nodel_pool:  candidates nodes for resection
% OUTPUTS:
%       new_idx:                indices of new target resection

    node_pool = data_model.node_pool;

    nr_nodes_in_list = numel(node_pool);
    nr_nodes_to_remove = numel(connection_indices);

    %%
    if nr_nodes_to_remove ~= nr_nodes_in_list

        % Select number of nodes to change from solution
        n_max_c = min(nr_nodes_to_remove, ...
            -nr_nodes_to_remove +nr_nodes_in_list);
        n_ch = powerlaw_range(1:n_max_c); % number of nodes to change

        % Idx of nodes to remove from VR list
        rm_idx = randperm(nr_nodes_to_remove, n_ch);
        % Idx of nodes to add to VR list
        add_idx = randperm(-nr_nodes_to_remove +nr_nodes_in_list, n_ch);

        %interchange index
        idx_not_VR = 1:nr_nodes_in_list;
        idx_not_VR(connection_indices) = [];

        % New connection indices
        new_idx = connection_indices;
        new_idx(rm_idx) = idx_not_VR(add_idx);

    else
        new_idx = connection_indices;
    end

    %check whether list consists of only unique numbers
    if numel(new_idx) ~= numel(unique(new_idx))
        error('Some connections are repeated in neighbourhood_fun\n');
    end