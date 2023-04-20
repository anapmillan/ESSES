function [order_in, order_out, frac_x] = func_SIR_adaptive_inverted(...
    aij, gamma, seed)
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
% function [order_in, order_out, frac_x] = ...
%           func_SIR_adaptive(aij, gamma)
% Simulate SIR spreading over a weighted, symmetric network for a given
% recovery probability gamma. The spreading probabilities are given by the
% network weights.
% This funciton uses the parallel Computing Toolbox. It is highly
% recommended to use to speed-up the SIR simulations.
%
% INPUTS:
%   aij:    adjacency matrix (weighted, symmetric)
%   gamma:  recovery probability
%   seed:   set of seed regions (either RA or derived from the
%           seed-probability maps)
% OUTPUTS:
%   order_in:  [nr_ROIs+1 x nr_ROIs] matrix where each element
%              (i,j) indicates the probability that ROI j becomes infected
%              at step i
%   order_out: [nr_ROIs+1 x nr_ROIs] matrix as order_in, but indicating the 
%              probability of recovery
%   frac_x:    [3 x nruns] matrix indicating the fraction of ROIs in the S,
%              I, and R state in the final step of the simulation.

    % Gather data:
    global data_model
    nr_ROIs = data_model.nrois;
    nruns   = data_model.nruns;
    order_in = zeros(nr_ROIs+1, nr_ROIs);
    %matrix to store order in which ROIs become active
    order_out = zeros(nr_ROIs+1, nr_ROIs);
    %matrix to store order in which ROIs become inactive

    % Output variable:
    frac_x = zeros(3, nruns);

    %% Prepare data
    % The code uses the list of neighbours of each node, stored in vvec, 
    % As well as the corresponding link weights, stored in vvec_w. 

    bij     = aij>0;        % Binary adjacency matrix
    ki      = sum(bij);     % Number of neighours of each node
    kmax    = max(ki);      % Maximum number of neighbours: used for vvec
    vvec    = zeros(nr_ROIs,kmax);   % Neighbours of each node
    vvec_w  = zeros(nr_ROIs,kmax);   % vvec_w(i,j) is the weight of the 
                                     % connection between node i and its 
                                     % j-th neighbour
    % Fill the neighbour matrices:
    for i = 1:nr_ROIs
        vvec(i,1:ki(i)) = find(bij(i,:));
        vvec_w(i,1:ki(i)) = aij(i,vvec(i,1:ki(i)));
    end

    %% Run SIR dynamics
    for ir = 1:nruns
        % Internal variables
        v_order     = zeros(nr_ROIs,1); % Activation order
        v_order_out = zeros(nr_ROIs,1); % Recovery order
        vvec_in     = vvec;
        vvec_w_in   = vvec_w;

        % Initial condition: RA is active
        vinf    = seed;         % vector of infected nodes
        v_order(vinf)   = 1;    
        vrec    = [];           % vector of recovered nodes
        t       = 2;            % counter: infection order

        % RECOVERY mechanism
        t_out   = 1;            % counter: recovery order
        % Infected nodes may recover with probability gamma:
        who_rec = rand(numel(vinf),1) < gamma;
        % If any nodes recover, update vectors:
        if (sum(who_rec)>0)
            % Recovery order:
            v_order_out(who_rec) = t_out;             
            t_out = t_out + 1;     
            % Vectors of recovered and infected nodes:
            vrec  = [vrec, vinf(who_rec)]; 
            vinf(who_rec)        = [];
        end
        % Vector of infected and recovered nodes
        vinft = [vinf, vrec]; 

        % INFECTION mechanism: TARGETS: neighbour of vinf that are not in
        % vinft
        % Targets:
        f1  = vvec_in(vinf,1:kmax)';    % identify neighbours of infected nodes
                                        % we count all, infection is per
                                        % connection with infected node
        f1  = f1(f1>0);
        % Weights:
        wf1 = vvec_w_in(vinf,1:kmax)';  % probability of infection of the neighbours
        wf1 = wf1(wf1>0);

        % Remove already inf/rec
        fr = f1(~ismember(f1,vinft));
        wfr = wf1(~ismember(f1,vinft));
        nfr = numel(fr);

        %DYNAMICS LOOP
        while (numel(vinf)>0 && nfr>0)  %while there are infected nodes and
                                        %the number of targets is >0
            %Pick new infected
            auxn    = sum(wfr) * rand;
            cwfr    = cumsum(wfr);
            nup     = min(find(cwfr>auxn));
            new_inf = fr(nup);
            % Store new infected in v_order and vinf, update counter
            v_order(new_inf) = t;
            t                = t+1;
            vinf             = [vinf, new_inf];

            % Recovery mechanism:
            who_rec = rand(numel(vinf),1) < gamma;
            if (sum(who_rec)>0)
                v_order_out(who_rec) = t_out;
                t_out = t_out + 1;
                vrec = [vrec, vinf(who_rec)];
                vinf(who_rec) = [];
            end
            vinft = [vinf, vrec];

            % Get TARGETs: neighrbous of infected nodes
            f1      = vvec_in(vinf,:);
            f1      = f1(f1>0);
            keep    = ~ismember(f1,vinft);
            fr      = f1(keep>0);
            nfr     = numel(fr);

            wf1     = vvec_w_in(vinf,:);
            wf1     = wf1(wf1>0);
            wfr     = wf1(keep>0);
        end


    %%
        % Store iteration results:
        order_in = order_in + ...
            sparse(v_order+1, 1:nr_ROIs, ones(1,nr_ROIs), 1+nr_ROIs,nr_ROIs);
        order_out = order_out + ...
            sparse(v_order_out+1, 1:nr_ROIs, ones(1,nr_ROIs), 1+nr_ROIs,nr_ROIs);

        frac_x(:,ir) = [(nr_ROIs -numel(vrec) -numel(vinf)),...
            numel(vinf), numel(vrec)];
    end

    % Normalize fraction of nodes:
    frac_x = frac_x/nr_ROIs;