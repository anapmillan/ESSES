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
% This script performs a virtual resection of the resection area. The SIR
% parameters are based on the best fit obtained with
% "spreading_model_fit.m" 
% The seed regions are based on the seed-probability map.

% This script makes use of the following patient-specific data:
%   1.- Brain network....................... file: data/network_example.txt 
%   2.- Seed regions (e.g. resection area).. file: data/seed_example.txt
%   3.- Seed probability map................ file: data/pattern_example.mat 
% The script outputs a figure with a violin plot of the spreading before
% and after the virtual resection.

%% Relevant data
nreps   = 100; % Number of iterations of the VR algorithm
nruns   = 1e4; % Number of SIR iterations in each repetition

% Load SIR parameters:
data_fit    = dlmread('data/best_fit.txt');
rhos        = data_fit(2); % Network density
gammas      = data_fit(3); % Recovery probability

tmax0   = 500;      % Maximum integration time (number of steps)
nrois   = 246;      % Number of nodes
seed_RS     = 3;    % Re-scaling parameters for the seed-probability maps
verbosity   = 0;

% Tags to get the patient data
get_SEEG    = 0;
get_map     = 1;

% Files for patient-specific inputs:
pat_dir = struct();
pat_dir.seed_filename       = 'data/seed_example.txt';      % File with seed data (read seed)
pat_dir.network_filename    = 'data/network_example.txt';   % File with network data (read wij)
pat_dir.pattern_filename    = 'data/pattern_example.mat';   % File with pattern data (load pattern)
pat_dir.map_filename        = 'data/seed_map_example.txt';  % File with seed probability map

%% Paths for outputs
pout_dyn        = 'VR_RA';          % SIR results
pout_figures    = 'VR_RA_results';  % VR results
% Create directories if they do not exist:
if ~isfolder(pout_dyn);     mkdir(pout_dyn);        end
if ~isfolder(pout_figures); mkdir(pout_figures);    end

%% Load paths and patient data
addpath('model')
addpath(genpath('toolbox'))

%Read network, resection area and seed probability map:
[wij, RA, ~, seed_map] = ...
    prepare_data_patient(pat_dir, get_SEEG, get_map);

% This tag will be used for the filenames:
name_tag = sprintf('inverted_seed_pw%.1f_%druns_%dreps', ...
    seed_RS, nruns, nreps);
% Name of final figure:
name_results = 'VR_RA_results';

%% Prepare variables
% The data for the simulations is stored in three structure variables: 
% data_model, data_patient and data_run
global data_model data_patient
data_model      = struct();
data_patient    = struct();
fill_data_variables; %save everything into data_model
data_model.w_seed       = 'map';

seed_map = seed_map.^seed_RS;
data_patient.nrois_soz  = sum(seed_map); %re-use variable
data_patient.esoz       = seed_map;
data_patient.RA         = RA;
data_patient.aij     = threshold_matrix(wij, data_run.km0);

data_run    = struct();
data_run.gamma  = gammas;
data_run.km0    = rhos*nrois;

%% Virtual resection of the RA:
fprintf('Running VR simulations...')
[spreading, seed_size] = run_VR(nreps, data_run);
save(sprintf('%s/spreading_%s', pout_dyn, name_tag), ...
    'spreading', 'seed_size');

%% Spreading stats
% We store the mean and std of 4 variables:
%   1.- Baseline spreading = IR(BS)
%   2.- Spreading after the VR = IR(VR)
%   3.- Decrease in spreading = IR(BS)-IR(VR)
%   4.- Normalized decrease in spreading = [IR(BS) - IR(VR)]/IR(BS)

[spreading_mean, spreading_std] = deal(zeros(1, 4));
spreading_mean(1:2) = mean(spreading);
spreading_std(1:2)  = std(spreading);

spreading_change    = spreading(1) - spreading(2);
spreading_mean(3)   = mean(spreading_change);
spreading_std(3)    = std(spreading_change);

spreading_rel_chage = spreading_change./spreading(:,1);
spreading_mean(4)   = mean(spreading_rel_chage);
spreading_std(4)    = std(spreading_rel_chage);

% Measure effect of resection:
spreading_diff      = [0 0];
spreading_diff(1)   = spreading_mean(3);
[~, spreading_diff(2)] = ttest(spreading(:,1), spreading(:,2));

fprintf('Baseline spreading IR(BS)..........= %.3f +- %.3f\n', ...
    spreading_mean(1), spreading_std(1))
fprintf('Spreading after resection IR(VR)...= %.3f +- %.3f\n', ...
    spreading_mean(2), spreading_std(2))
fprintf('Decrease in spreading..............= %.3f +- %.3f (p=%.2f)\n', ...
    -spreading_mean(3), spreading_std(3), spreading_diff(2))
fprintf('Relative Decrease in spreading.....= %.3f +- %.3f\n', ...
    -spreading_mean(4), spreading_std(4))

%% Reformat data for violin plot
spreading_all   = zeros(nreps*2,1);  % IR_BS, IR_VR (300 reps)
group_all       = strings(nreps*2, 1);
str_VR          = {'BS','VR'};
for j   = 1:2
    i0  = (j-1)*nreps +1;
    i1  = j*nreps;
    spreading_all(i0:i1)    = spreading(:,j);
    group_all(i0:i1)        = repmat(str_VR{j}, nreps, 1);
end
group_all   = cellstr(group_all);
 
%% Violin plot of spreading before and after resection
figure
ccc = [[1 0 0]; 0.7*[1 0 0]];
vs  = violinplot(spreading_all, group_all, 'ViolinColor',ccc);
for i   = 1:numel(vs)
    vs(i).ScatterPlot.SizeData  = 28;
    vs(i).ScatterPlot2.SizeData = 28;
end
ylabel('$IR$', 'Interpreter', 'latex')  
set(gca, 'FontSize', 14, 'FontName', 'Times')

set(gcf, 'paperunits','inches','paperpositionmode','manual',...
'papersize',[5 3], 'paperposition', [0 0 5 3])
print('-dpng', sprintf('%s/%s.png', pout_figures, name_results))
 