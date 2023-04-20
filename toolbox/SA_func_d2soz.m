function [rmin_sol, set_rm_SA, VR_sizes] = ...
    SA_func_d2soz(w_use, rec_set, cost1)
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
%   function [rmin_sol, set_rm_SA, VR_sizes] = ...
%       SA_func_d2soz(w_use, rec_set, cost1)
%   Implement SA algorithm to identify optimal resections of increasing
%   size.
%   INPUTS:
%       w_use:      String to select initial condition. If 'PREV', uses the
%           S-1 solution plus a randomly selected node as the initial 
%           condition to find the alternative resection of size S.
%       rec_set:    Initial condition for S=1
%       cost1:      Cost for S=1
%   IMPLICIT INPUTS:
%       data_model,     with field node_pool
%       data_patient,   with field aij
%       data_SA
%   OUTPUTS:
%       rmin_sol:   Efficiency of alternative resections
%       set_rm_SA:  Alternative resections
%       VR_sizes:   Sizes of alternative resections

    fprintf('Running SA\n') 
    %% Initialisation
    global data_model data_SA data_patient

    % Load data
    aij             = data_patient.aij;                 % Original adj. matrix
    repeat          = data_SA.repeat;                   % No. SA iterations
    % Possible target nodes
    node_pool       = data_model.node_pool;             
    n_nodes_pool    = numel(node_pool); 
    % VR sizes
    VR_sizes        = data_SA.VR_sizes;                  
    n_VR_sizes      = numel(VR_sizes);
    % Seed regions:
    esoz            = find(data_patient.esoz);
    % Seed-probability map
    esoz_map        = data_patient.esoz;
    % Auxiliary variables for the parfor loop:
    data_SA_aux     = data_SA;
    data_model_aux  = data_model;

    % Define output variables:
    rmin_sol    = zeros(size(VR_sizes,2),repeat);
    cost        = zeros(n_VR_sizes, repeat);
    set_rm_SA   = cell(repeat, 1);

    % Function handles
    handle_neighbourhood    = @(x) neighbourhood_func(x, data_model_aux);
    handle_loss             = @(x) loss_func(x, aij, esoz_map);
 

    %% Run SA
    tic
    parfor r   = 1:repeat   % Perform n = repeat iterations of the SA optimization
                            % recursively from size S=1 until the
                            % distonnection size
        % Define internatl variables and initial condition                    
        node_pool_aux   = node_pool;
        x0_aux          = rec_set;
        w_use_aux       = w_use;
        set_rm_SA_r     = cell(n_VR_sizes,1);
        VR_sizes_aux    = VR_sizes;

        set_rm_SA_r{1}  = rec_set(1);
        cost_in         = zeros(n_VR_sizes,1);
        cost_in(1)      = cost1;

        % Loop over resection sizes:
        for n = 2:n_VR_sizes
            nr_rem_connections = VR_sizes_aux(n);

            if  strcmp(w_use_aux, 'PREV')   % Use previous solution as initial condition
                set0        = set_rm_SA_r{n-1};
                rem_set     = 1:246;
                rem_set(set0)   = [];
                imax    = 0;
                if r<6 % Use best 1-node option for 5 runs then test random
                    d2soz_in    = zeros(numel(rem_set),1);
                    for i = 1:(246-n+1)
                       seti = [set0 rem_set(i)];
                       d2soz_in(i) = VR_d2soz(seti, aij, esoz);
                    end
                    [~, imax] = max(d2soz_in);
                else
                    imax    = randi(numel(imax));
                end
                initial_connection_indices = [set0 rem_set(imax)];
            else % Random initial condition
                initial_connection_indices = ...
                    randperm(n_nodes_pool, nr_rem_connections);
            end
            
            % Run simulated annealing:
            [final_indices, cost_in(n)] = simulated_annealing(...
                handle_neighbourhood, handle_loss,...
                initial_connection_indices, data_SA_aux);
            list_rm = node_pool_aux(final_indices);
            set_rm_SA_r{n} = list_rm;

            if cost_in(n) == 0
                fprintf('Seed disconnected for n=%d, r=%d\n', n, r)
                break;
            end
        end
        set_rm_SA{r} = set_rm_SA_r;
        cost(:,r)   = cost_in;
    end
    toc

    rmin_sol(1:size(VR_sizes,2),:) = cost;

end