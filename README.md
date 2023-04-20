# ESSES: Epidemic Spreading Seizure and Epilepsy Sugery framework.

This repository contains the MATLAB toolbox for the ESSES framework as defined in [1]. ESSES creates patient-specific seizure propagation models to investigate alternative resection strategies and estimate the outcome of a given resection for patients canditates for epilepsy surgery.
Seizure propagation is simulated by means of the SIR (Susceptible-Infected-Recovered) model, run on top of patient-specific brain connectivity. 
An early version of this toolbox can be found in https://github.com/anapmillan/computer_model_for_epilepsy

This toolbox allows for three analyses:
1.- Model definition and fit to data (spreading_model_fit.m):
Simulates seizure propagation for a set of parameters of the SIR model. Fits the model parameters to seizure propagation data derived from invasive EEG recordings. 
Produces a fit file ("data/best_fit.txt") with the parameters values (network density and recovery probability) leading to the best fit.

2.- Simulation of the effect of a given resection ("virtual_resection_RA.m").
Simulates the effect of a given resection strategy on seizure propagation.

3.- Optimization of alternative resection strategies.
Searches for patient-specific alternative resection strategies using simmulated annealing.


 The toolbox makes use of four inputs included under "data":
 
 1.- Network structure (network_example.txt)
 
     Backbone for seizure propagation. An artifitial exemplary network based on the Exponential Distance Rule on a 3D random distribution of centroids (to simulate a brain network) is included. The size has been set to N=246 in accordance to the Braintome Atlas. The script used to create this network can be found in https://github.com/anapmillan/computer_model_for_epilepsy
     
 2.- Epidemic seed (seed_example.txt)
 
     Seed or origin for the epidemics in the model. An artificial exemplary seed is included.
     
 3.- Seizure pattern (pattern_example.mat)
 
     Spreading pattern indicating the activation order of the different brain regions during a seizure. This a structure "patttern" including the fields:
     
     pattern.trois:      array indicating the indices of sampled active ROIs
     pattern.all_erois:  (n_sampled_rois x 1) array indicating all sampled ROIs
     pattern.order:      array indicating the activation order of active sampled
                         ROIs. Size = (n_sampled_active_rois x 1)
     pattern.ntrois:     total number of sampled ROIs
     
 4.- Seed-probability map (seed_map_example.txt).
 
    Probability of each ROI of starting a region. Can be derived from the presurgical evaluation of the patient. An artifitial exemplary map is included.
    
    

These programs are distributed by the authors in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

If you use any of these codes, then please cite:
[1] Ana P. Millan, et al. "Individualized epidemic spreading models predict epilepsy surgery outcomes: a pseudo-prospective study."  medRxiv (2023): 2023-03. doi: https://doi.org/10.1101/2023.03.16.23287370
[2] Ana P. Millan, et al. "The role of epidemic spreading in seizure  dynamics and epilepsy surgery." Network Neuroscience (2022): 1-55. doi: https://doi.org/10.1162/netn_a_00305

The toolbox makes use of the anneal function by Kirkpatrick et al. If you use it, please cite:
Kirkpatrick, S., Gelatt, C.D., & Vecchi, M.P. (1983). Optimization by Simulated Annealing. Science, 220, 671-680.
Contact: joachim.vandekerckhove@psy.kuleuven.be

To generate the figures, the following toolbox is needed as well:
https://github.com/bastibe/Violinplot-Matlab, DOI: 10.5281/zenodo.4559847



(c) Ana P. Millan (a.p.millanvidal@amsterdamumc.nl)

