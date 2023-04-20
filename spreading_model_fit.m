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

%% 
% This script simulates the SIR model for a set of parameters to find the 
% parameter configuration that best approximates iEEG-recorded seizure
% propagation. 

% The SIR parameters are the network density (rho) and the recovery
% probability (gamma). The spreading probability is set directly by the
% network density. 
% The seed regions are given by the resection area, used as a proxy for the
% seizure-onset-zone.

% This script makes use of the following patient-specific data:
%   1.- Brain network....................... file: data/network_example.txt 
%   2.- Seed regions (e.g. resection area).. file: data/seed_example.txt
%   3.- iEEG seizure propagation pattern.... file: data/pattern_example.mat 
% The script outputs a figure with the result fit map and the coordinates 
% (rho, gamma, fit_value) of the best fit


%% Relevant data
% Define mesh of gamma and density (rho) values
% Both are in log-scale
gammas  = [0.005, logspace(-2,0,10)];
rhos    = [logspace(-2,-0.7,6), 0.35];

nreps       = 10;    % number of repetitions to average C
nruns       = 1e4;  % number of SIR iterations in each repetition
get_SEEG    = 1;    % Flag to retrieve SEEG data
get_map     = 0;
verbosity   = 0;    % Verbosity level for output

% Files for patient-specific inputs:
pat_dir = struct();
pat_dir.seed_filename       = 'data/seed_example.txt'; % File with seed data (read seed)
pat_dir.network_filename    = 'data/network_example.txt'; % File with network data (read wij)
pat_dir.pattern_filename    = 'data/pattern_example.mat'; % File with pattern data (load pattern)
pat_dir.map_filename        = 'data/seed_map_example.txt'; % File with seed probability map

%% Paths for outputs
pout_dyn        = 'spreading_data';     % SIR data
pout_corr       = 'spreading_output';   % Model fit results
pout_figures    = 'model_fit_results';  % Results

% Create directories if they do not exist:
if ~isfolder(pout_dyn);     mkdir(pout_dyn);        end
if ~isfolder(pout_corr);    mkdir(pout_corr);       end
if ~isfolder(pout_figures); mkdir(pout_figures);    end

%% Parameters for model
tmax0   = 500;           % Maximum number of integration steps
nrois   = 246;           % Number of ROIs
kmos    = rhos*nrois;    % Mean degree
nks     = numel(kmos);   % Number of different densities
ngammas = numel(gammas); % Number of gamma values

%% Load paths and patient data
% Add necessary paths to the model and toolbox directories
addpath('model')
addpath(genpath('toolbox'))

% Prepare network data, resection area (seed), and iEEG pattern
[wij, RA, pattern] = prepare_data_patient(pat_dir, get_SEEG, get_map); 
seed = RA';

% Define name for correlation results file and figure
name_results = sprintf(...
    '%s/wcorr_%dgammas_%drhos_%druns',...
    pout_figures, numel(gammas), numel(rhos), nruns);

%% Prepare variables 
% The data for the simulations is stored in three structure variables: 
% data_model, data_patient and data_run

global data_model data_patient  
data_model      = struct();
data_patient    = struct();
fill_data_variables; %save everything into data_model
data_model.pout_corr    = pout_corr;
data_model.w_seed       = 'RA';
data_model.map          = 0;
data_patient.nrois_soz = numel(RA);
data_patient.esoz = RA;
data_patient.pattern = pattern;

%% Simulations
fprintf('Running SIR simulations...')
for itest = 1:nreps % LOOP over repetitions
    data_run        = struct();
    data_run.rep    = itest;
    data_run.is_VR  = '';
    data_patient.name_corr_file = sprintf(...
        '%s/wcorr_%s_%dgammas_%drhos_%druns_%d.mat',...
        pout_corr, numel(data_model.gammas),...
        numel(data_model.rhos), data_model.nruns, itest);

    fprintf('\n REP (%d/%d)\n', itest, nreps); 
    tic    
    for ik  = 1:nks % LOOP over network densities
        km0 = kmos(ik);
        data_run.km0 = km0;
        if verbosity;s0 = '\n'; else; s0 = ''; end
        fprintf('%s\t km = %.2f (%d/%d)\n',s0,km0, ik, nks);

        % Threshold matrix at required density
        aij = threshold_matrix(wij,km0);
        
        for igamma = 1:numel(gammas) % LOOP over recovery values
            if verbosity
                fprintf('g = %.2f (%d/%d)\n',...
                    gammas(igamma), igamma, numel(gammas));
            end
            data_run.gamma  = gammas(igamma);
            tot_pov         = spreading_func(aij, data_run, seed);
        end
        if verbosity; toc; end
    end
    toc
%%  Measure model fit
    output = correlation_analysis(itest); 
    if itest==1
        fit_map = output.pcp_weighted.*...
            output.frac_act_w.*output.frac_inact_w;
    else
        fit_map = fit_map +...
            output.pcp_weighted.*...
            output.frac_act_w.*output.frac_inact_w;
    end
end

% Save fit results
fit_map = fit_map/nreps;
save(sprintf('%s/%s', pout_figures, name_results), 'fit_map');

%% 
% This figure plots the fit results.
figure
imagesc(fit_map)  
set(gca, 'FontSize', 14, 'FontName', 'Times')
hold on
[i1,i2] = max(fit_map(:)); %best fit
[i2,i3] = ind2sub(size(fit_map), i2);
best_fit_ind = [i1 i2 i3];
plot(i3,i2, '*r', 'LineWidth',2)
colorbar
caxis([0 i1])
set(gca, 'XTick',1:numel(rhos), 'XTickLabel', num2str(rhos', '%.3f'))
set(gca, 'YTick',1:2:numel(gammas), ...
    'YTickLabel', num2str(gammas(1:2:end)', '%.2f'))
xlabel('Network density')
ylabel('Recovery probability')

set(gcf, 'paperunits','inches','paperpositionmode','manual',...
'papersize',[5 3], 'paperposition', [0 0 5 3])
print('-dpng', sprintf('%s/%s.png', pout_figures, name_results))

%% Save data of best fit
best_rho    = rhos(i3);
best_gamma  = gammas(i2);
best_fit    = i1;
dlmwrite('data/best_fit.txt', [best_fit, best_rho, best_gamma])

