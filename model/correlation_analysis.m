function output = correlation_analysis(irep)
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
% function output = correlation_analysis(irep)
%   Measure the goodness-of-fit of results from iteration irep.
%   The target filename must match the one in spreading_model_fit.m
% INPUT: 
%   irep: Iteration number. Used for the filenames.
% IMPLICIT INPUTS:
%   Requires the following fields in 
%   data_model:   nrois, pout_dyn, pout_corr, gammas, rhos, nruns
%   data_patient: pattern
% OUTPUTS:
%   output: Structure variable with the fields:
%       gammas:       [1 x n_gammas] matrix equal to data_model.gammas;
%       rhos:         [1 x n_rhos] matrix equal data_model.rhos;
%       pcp_weighted, frac_eq_w, n_act_model, frac_act_w, frac_inact_w:
%                     [n_gammas x n_rhos] matrices indicating, for each 
%                     (gamma,rho) point, the:
%                       - weighted Pearson Correlation Coefficient 
%                       - weighted fraction of ROIs in the same state (active/inactive)
%                       - weighted fraction of ROIs that are active 
%                       - weighted fraction of ROIs that are inactive 
%                     between iEEG and SIR seizure propagation patterns for each 
%                     (gamma,rho) point.
%       n_act_model:  [n_gammas x n_rhos] indicating the number of active ROIs 
%                     in the SIR patterns, for each (gamma,rho) point
   

    global data_model data_patient

    % Gather data
    pattern     = data_patient.pattern;
    nrois       = data_model.nrois;
    pin         = data_model.pout_dyn;
    n_gammas    = numel(data_model.gammas);
    n_rhos      = numel(data_model.rhos);

    % Create output variables:
    [frac_eq_w, n_act_model, pcp_weighted, frac_act_w, frac_inact_w] = ...
        deal(zeros(n_gammas, n_rhos));

    %% Identify Sampled, active, inactive ROIs
    % All sampled ROIs:
    all_samp_rois   = unique(pattern.all_erois);
    all_samp_rois   = all_samp_rois(all_samp_rois>0);
    nrois_sampled   = numel(all_samp_rois);

    % Sampled and active (iEEG) ROIs:
    r_rois          = pattern.trois;
    %Remove any repeated entries:
    [aux,iaux]  = unique(r_rois);
    r_rois      = r_rois(sort(iaux));

    % Sampled and inactive (iEEG) ROIs:
    r2_rois     = setdiff(all_samp_rois, r_rois); %Inactive ROIs in pattern
    r2_rois     = r2_rois(r2_rois>0); 

    %% Loops over files
    for igamma = 1:n_gammas
        for irho = 1:n_rhos
            %% Read SIR data
            name_data = sprintf('%s/dyn___km%.4f_gamma%.4f_nruns%d_rep%d.txt',...
                pin, 246*data_model.rhos(irho),...
                data_model.gammas(igamma), data_model.nruns, irep);
            dd  = dlmread(name_data);   % read file
            dd  = dd/max(dd(:));        % normalize

            % Measure mean activation time from distribution
            nd2     = size(dd,2);
            auxn    = sum(dd(:,2:end),2);
            auxn(auxn==0) = 1;
            dd0     = dd(:,1);
            dd      = dd./auxn;
            dm      = dd*[0:(nd2-1)]'; %mean activation time

            %% Measure overlap between the iEEG and SIR patterns

            % Active ROIs in model (m_rois) and both (mr_rois)
            m_rois  = find(dm);
            mr_rois = intersect(m_rois,r_rois); 

            % Mesure overlap:
            frac_eq_w(igamma, irho)     = (sum(1-dd0(r_rois)) +...
                sum(dd0(r2_rois))) / nrois_sampled;
            % (Probability of being active of active rois + prob. of being
            % inactive of inactive rois ) / number of sampled rois
            frac_act_w(igamma, irho)    = mean(1-dd0(r_rois));
            frac_inact_w(igamma, irho)  = mean(dd0(r2_rois));

            % Number of active ROIs in SIR pattern:
            n_act_model(igamma, irho)   = numel(m_rois)/nrois;

            % Weights for correlation
            wij_roi     = ones(numel(mr_rois,1),1) - dd0(mr_rois); %active


            %% Measure weighted PCP correlation bewteen SIR and iEEG patterns
            % The correlation is measured over the activation orders of the 
            % sampled and active ROIs (in both patterns):

            % Activation in SIR model
            dm  = dm(mr_rois); %keep only sampled and active rois
            dmi = unique(dm);
            i0  = 0;
            act_order_SIR  = zeros(size(dm));
            for iid = 1:numel(dmi)
                idx = dm==dmi(iid);
                act_order_SIR(idx) = i0+mean(1:sum(idx));
                i0  = i0 + sum(idx);
            end

            %Activation in iEEG data
            datam = pattern.order; 
            ikeep = ismember(r_rois,mr_rois); %keep only active and sampled
            datam = datam(ikeep);
            act_order_iEEG = zeros(size(datam));
            datamu = unique(datam);
            i0 = 0;
            for iid = 1:numel(datamu)
                idx = datam==datamu(iid);
                act_order_iEEG(idx) = i0 + mean(1:sum(idx));
                i0 = i0 + sum(idx);
            end

            % Measure correlation only if there are enough points
            if (numel(act_order_SIR) <2 | ...
                numel(unique(act_order_SIR))<2 | ...
                numel(unique(act_order_iEEG))<2 )
                pcp_weighted(igamma, irho) = 0;
            else
                aux = weightedcorrs([act_order_iEEG,act_order_SIR], wij_roi');
                pcp_weighted(igamma, irho) = aux(2);
            end


        end
    end

    %% Create and fill output variable
    %It is used to store the results and returned to the main progra
    output = struct();
    output.gammas   = data_model.gammas;
    output.rhos     = data_model.rhos;
    output.pcp_weighted = pcp_weighted;
    output.frac_eq_w    = frac_eq_w;
    output.n_act_model  = n_act_model;
    output.frac_act_w   = frac_act_w;
    output.frac_inact_w = frac_inact_w;

end