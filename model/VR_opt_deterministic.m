function [d2soz0, d2soz_lin, rm_edges_sorted, d2soz_rec, rec_set] = ...
    VR_opt_deterministic()
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
% Perform deterministic search of alternative resection strategies.
% IMPLICIT INPUTS:
%   data_model.node_pool:   target nodes for VRs
%   data_patient.aij:       [Nr_ROIs x Nr_ROIs] unresected adjacency matrix
%   data_patient.esoz:      [Nr_ROIs x 1] seed-probability map
% OUTPUTS:
%   d2soz0:             Seed efficiency in the orignal network 
%   d2soz_lin:          [Nr_ROIs x 1] vector indicating the seed efficiency
%                       following 1-node VRs.
%   rm_edges_sorted:    Indices of nodes sorted according to (minimum)
%                       seed efficiency following their VR.
%   d2soz_rec:          [Nr_ROIs x 1] vector indicating the seed efficiency
%                       following the best recursive resection of size j.
%                       It is 0 by definition for all sizes greater than
%                       the disconnecting size.
%   rec_set:            Indices of nodes that are added to the recursive
%                       resection at each size, up to the size of the 
%                       disconnecting resection. Thus, nonzeros(rec_set) is 
%                       the disconnecting resection. 

    global data_model data_patient
    aij         = data_patient.aij;
    node_pool   = data_model.node_pool;
    max_size_VR = numel(node_pool);
    esoz_map    = data_patient.esoz;

%% Base distance
    d2soz0 = matrix2distance(aij, esoz_map);

%% NOW LINEAR
    d2soz_lin = zeros(max_size_VR,1);
    fprintf('Running linear\n')
    tic
    parfor ir = 1:max_size_VR
        d2soz_lin(ir) = VR_d2soz(node_pool(ir), aij, esoz_map);
    end
    toc

%% Recursive deterministic
    fprintf('Running recursive\n')
    tic
    [d2soz_lin_sorted, rm_edges_sorted] = sort(d2soz_lin,'descend');
    rec_set     = zeros(max_size_VR,1);
    rec_set(1)  = rm_edges_sorted(1); %first edge: min ind corr.

    d2soz_rec       = zeros(max_size_VR,1);
    d2soz_rec(1)    = d2soz_lin_sorted(1);
    rem_set         = rm_edges_sorted(2:end);

    i       = 2;
    soz_con = 1;

    while i<=max_size_VR & soz_con
        % Test all possible links
        aux_d2soz = zeros(numel(rem_set),1);
        parfor ir = 1:numel(rem_set)
            node_pool_aux   = node_pool;
            rec_set_aux     = rec_set;
            rec_set_aux     = rec_set_aux(1:i-1);

            new_set_VR      = node_pool_aux([rec_set_aux; rem_set(ir)]);
            aux_d2soz(ir)   = VR_d2soz(new_set_VR, aij, esoz_map);
        end

        if any(aux_d2soz==0)
            fprintf('SOZ completely disconnected at link %d\n',i);
            aux0            = find(aux_d2soz==0);
            rec_set(i)      = rem_set(aux0(1));
            d2soz_rec(i)    = 0;
            soz_con         = 0;
        else
            [c_min, e_min]  = min(aux_d2soz);
            d2soz_rec(i)    = c_min;
            rec_set(i)      = rem_set(e_min);
            rem_set(e_min)  = [];
            i = i+1;
        end
    end
    toc