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
% This script optimizes virtual resections of increasing sizes. To optimize
% the resection strategy, the script uses the seed efficiency as a
% surrogate for SIR spreading. The effect of each resection is given by
% e(VR) = seed_efficiency(BS) - seed_efficiency(VR).
% iterative procedure is followed: 
%   1.- Linear aproximation: 1-node resections are performed and the nodes
%                            are ordered according to e(VR).
%   2.- Recursive aproximation: Starting form the most efective 1-node
%                               resection, the size of the resection is 
%                               increased recursively by adding the node 
%                               that maximally decrease seed_efficiency(VR).
%   3.- Simulated annealing: Identify alternative resections with using
%                            simmulated annealing (SA). The SA algorithm is 
%                            iterated nreps_SA times and the best solution 
%                            selected.
%   4.- Actual effect (SIR): The actual effect of each alternative
%                            resection is measured with the SIR model. The 
%                            resection leading to the maximum decrease in 
%                            spreading is selected, for each size.  
% For the definition of the seed regions in the surrogate analyses, weak
% seeds (lower probability than maximum -2*standard deviation of the seed)
% probability map) are included.

% The SIR parameters are based on the best fit obtained with
% "spreading_model_fit.m" 
% The seed regions are based on the seed-probability map.

% This script makes use of the following patient-specific data:
%   1.- Brain network....................... file: data/network_example.txt 
%   2.- Seed regions (e.g. resection area).. file: data/seed_example.txt
%   3.- Seed probability map................ file: data/pattern_example.mat 
% The script outputs a figure with a violin plot of the spreading before
% and after the virtual resection.


%% Relevant data
nreps_seed  = 100; % Number of seed realizations for each analysis.
nruns       = 1e4; % Number of SIR iterations in each repetition.
nreps_SA    = 40;  % Number of alternative SA resections.

network_scale   = 'log'; % Use log-scaling to measure seed efficiency.              
    
% Load SIR parameters:                              
data_fit    = dlmread('data/best_fit.txt');
rhos        = data_fit(2);
gammas      = data_fit(3);

tmax0   = 500;     % Maximum integration time (number of steps).
nrois   = 246;     % Number of nodes.
w_RS    = 3;       % Re-scaling parameter for the seed-probability maps.
verbosity   = 0;

% Tags to get patient data:
get_SEEG    = 0;
get_map     = 1;

% Files for patient-specific inputs:
pat_dir = struct();
pat_dir.seed_filename       = 'data/seed_example.txt';     % File with seed data (read seed)
pat_dir.network_filename    = 'data/network_example.txt';  % File with network data (read wij)
pat_dir.pattern_filename    = 'data/pattern_example.mat';  % File with pattern data (load pattern)
pat_dir.map_filename        = 'data/seed_map_example.txt'; % File with seed probability map


%% Paths for outputs
pout_dyn        = 'op_VR';         % SIR data
pout_figures    = 'op_VR_results'; % Results
% Create directories if they do not exist:
if ~isfolder(pout_dyn);      mkdir(pout_dyn);         end
if ~isfolder(pout_figures);  mkdir(pout_figures);    end

%% Load paths and patient data
addpath('model')
addpath(genpath('toolbox'))

%Read network, resection area and seed probability map:
[wij, RA, ~, seed_map] = ...
    prepare_data_patient(pat_dir, get_SEEG, get_map);

% This tags will be used for the filenames:
name_tag    = sprintf('VRop_seed_pw%d_%s_%dSAreps',...
    w_RS, network_scale, nreps_SA);
name_tag_SIR    = sprintf('%s_%dnreps_seed_%dnrunsSIR_test', ...
    name_tag, nreps_seed, nruns);


%% Prepare variables
% The data for the simulations is stored in four structure variables: 
% data_model, data_patient, data_SA and data_run

global data_model data_patient data_SA
data_model      = struct();
data_patient    = struct();
fill_data_variables; %save everything into data_model
data_model.w_seed       = 'map';
data_model.node_pool    = 1:246; 

% For surrogate metric, remove very weak seeds from seed pool
esoz    = seed_map.^seed_RS;
aux     = esoz;
aux(aux< (max(aux)-2*std(aux)))   = 0;
data_patient.nrois_soz  = sum(aux); %re-use variable
data_patient.esoz       = aux;
nrois_soz               = numel(aux);

% Get network
aij0    = threshold_matrix(wij, data_run.km0);
% For the surrogate metric we create a distance matrix instead
% Since the seed efficiency is used as a proxy
aij          = 1./aij0; % for the surrogate metric
aij(aij0==0) = 0; 
% Use log-scale for distance matrix:
if strcmp(network_scale, 'log')
    amin    = min(aij(aij>0));
    wijaux  = log10(aij/amin);
    wijaux(aij==0) = 0;
    aij     = wijaux/max(wijaux(:));
end
data_patient.aij     = aij; 

% Data for simulated annealing:
data_SA         = struct();
data_SA.repeat  = nreps_SA;
data_SA.VR_sizes    = 1:246;     
data_SA.w_netw  = network_scale; %raw or log
data_SA.options.InitTemp    = 1;
data_SA.options.StopTemp    = 1e-8;          
data_SA.options.Verbosity   = 0;
data_SA.options.CoolSched   = @(T) (.8*T);    
data_SA.options.MaxConsRej  = 1000;          

% Data for simulations
data_run        = struct();
data_run.gamma  = gammas;
data_run.km0    = rhos*246;

%% Output files
nout_VRop_det   = sprintf(...
    '%s/VRop_det_results_%s_%dsizes_%dto%d_%dreps.mat',...
    pout_dyn, name_tag, numel(data_SA.VR_sizes),...
    data_SA.VR_sizes(1), data_SA.VR_sizes(end), data_SA.repeat);
nout_VRop_SA    = sprintf(...
    '%s/VRop_SA_results_%s_%dsizes_%dto%d_%dreps.mat',...
    pout_dyn, name_tag, numel(data_SA.VR_sizes),...
    data_SA.VR_sizes(1), data_SA.VR_sizes(end), data_SA.repeat);

%% Surrogate analyses
% Deterministic analysis
[d2soz0, d2soz_lin, rm_nodes_sorted, d2soz_rec, rec_set] = ...
    VR_opt_deterministic();
save(nout_VRop_det, 'd2soz0','d2soz_lin','rm_nodes_sorted',...
    'd2soz_rec', 'rec_set')

% SA
[d2soz_SA, rm_nodes_SA, VR_sizes] = SA_func_d2soz('PREV',...
    rec_set, d2soz_rec(1));
save(nout_VRop_SA, 'd2soz_SA','rm_nodes_SA')
d2soz_SA    = d2soz_SA';


%% Plot results
% Plot effect of each alternative resection, with each method
ylabels = {'Individual Cont.', 'Recursive', 'SA'};
cc = parula(3);
figure 
hold on 
xlim([1 246])
xlabel('ROIs')
ylabel('Seed Efficiency')
 
plot(sort(d2soz_lin), ...
    'color',cc(1,:), 'linewidth',2)
plot(d2soz_rec, ...
    'color',cc(2,:), 'linewidth',2)
plot(data_SA.VR_sizes, squeeze(d2soz_SA)', ...
    'color',cc(3,:), 'linewidth',0.5, 'HandleVisibility','off')
plot(data_SA.VR_sizes, squeeze(min(d2soz_SA, [], 1)), ...
    'color',cc(3,:), 'linewidth',2)
legend(ylabels, 'Location','southeast')

set(gca, 'FontSize', 14, 'FontName', 'Times')
set(gca,'XScale', 'log', 'YScale', 'log')
set(gcf, 'paperunits','inches','paperpositionmode','manual',...
'papersize',[5 3], 'paperposition', [0 0 5 3])
print('-dpng', sprintf('%s/optimization_analysis_%s', ...
    pout_figures, name_tag))
 

%% Measure actual effect of the VR (best SA only)
% Update data for plots and simulations:
nout_SA_SIR = sprintf('%s/run_VRop_data_%s.mat', ...
    pout_dyn, name_tag_SIR);
data_model.nruns    = nruns;
data_model.name_tag = name_tag_SIR;
data_patient.aij = aij0;
data_patient.esoz = esoz;
data_patient.nrois_soz = numel(esoz);

% Output variables:
spreading_VR = zeros(nreps_seed, numel(data_SA.VR_sizes));
spreading_BS = zeros(nreps_seed, 1);
seed_size    = zeros(nreps_seed,1);

% Best SA iteration:
[aux, idx]  = min(squeeze(d2soz_SA));
max_size    = min(find(aux==0));

fprintf('Running SIR dynamics...')
for irep    = 1:nreps_seed
    fprintf('%d ', irep)
    data_run.rep    = irep;

    %Get new seed: 
    while seed_size(irep)==0
        new_seed    = find(data_patient.esoz>rand(246,1))';
        seed_size(irep)   = numel(new_seed);
    end 

    % Baseline spreading:
    data_run.is_VR      = 'BS'; %tag for the namefiles so that BS and VR are diff
    spreading_BS(irep)    = mean_SIR_spread(aij, data_run, new_seed);

    % Perform VRs of increasing sizes:
    for j = 1:max_size
        aij_VR = prune_nodes(aij, rm_nodes_SA{idx(j)}{j});
        data_run.is_VR  = sprintf('VR%d', data_SA.VR_sizes(j));
        %tag for the namefiles so that BS and VR are diff
        % and to differentiate VRs of diff sizes
         spreading_VR(irep, j) = ...
             mean_SIR_spread(aij_VR, data_run, new_seed);
    end
end
save(sprintf('%s/SIR_results.mat', pout_dyn),...
    'spreading_VR', 'spreading_BS', 'seed_size')

%% Plot SIR results: all iterations
cc = parula(1);
cc_i = 0.5*[1 1 1 1];
figure
hold on
for irep = 1:nreps_seed
    plot(data_SA.VR_sizes, squeeze(spreading_VR(irep, :)), '.', ...
        'color', cc_i, 'MarkerSize', 2);
end
xlim([1 max_size])
ylabel('Spreading')
xlabel('VR Size')
plot(data_SA.VR_sizes, squeeze(mean(spreading_VR,1)), '-', ...
    'linewidth',2, 'color', cc);

set(gca, 'FontSize', 14, 'FontName', 'Times')
set(gcf, 'paperunits','inches','paperpositionmode','manual',...
'papersize',[5 3], 'paperposition', [0 0 5 3]) 
print('-dpng', sprintf('%s/VR_effect_points_AV_%s.png', ...
    pout_figures, name_tag_SIR));

%% Measure average and relative change
mean_spreading_VR   = squeeze(mean(spreading_VR, 1));
mean_spreading_BS   = mean(spreading_BS);
std_spreading_VR    = squeeze(std(spreading_VR, [], 1));
dec_spreading_VR    = mean_spreading_BS - mean_spreading_VR;
ndec_spreading_VR   = dec_spreading_VR./mean_spreading_BS;

%% Plot SIR results: average and decrease
ylabels = {'$IR_{BS}$', '$IR_{VR}$', '$\Delta IR$', '$\delta IR$'};
cc = parula(3);
figure 
hh1 = subplot(3,1,[1 2]);
hold on
box on
xlim([1 max_size])
xlabel('$S_{VR}$', 'interpreter','latex') 
title('Effect in spreading')
plot(data_SA.VR_sizes, mean_spreading_BS*ones(size(data_SA.VR_sizes)),...
    '--k', 'linewidth',2)
plot(data_SA.VR_sizes, mean_spreading_VR,...
    'color', cc(1,:), 'linewidth',2)
plot(data_SA.VR_sizes, dec_spreading_VR,...
    'color', cc(2,:), 'linewidth',2)
set(gca, 'FontSize', 14, 'FontName', 'Times')
legend(ylabels(1:3), 'Interpreter','latex', 'location','northwest')

hh2 = subplot(313);
hold on
box on
xlim([1 max_size])
xlabel('$S_{VR}$', 'interpreter','latex') 
plot(data_SA.VR_sizes, ndec_spreading_VR,...
    'color', cc(3,:), 'linewidth',2)
set(gca, 'FontSize', 14, 'FontName', 'Times')
legend(ylabels(end), 'Interpreter','latex', 'location','northwest')


set(gcf, 'paperunits','inches','paperpositionmode','manual',...
'papersize',[5 3], 'paperposition', [0 0 5 3])  
print('-dpng', sprintf('%s/VRop_effect_av_curves_%s.png', ...
    pout_figures, name_tag_SIR));