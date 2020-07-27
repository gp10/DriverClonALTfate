%% SCRIPT USED FOR ANALYZING WT & MUTANT TOTAL CLONE DYNAMICS:
% WT and Mutant clone dynamics are simulated using a Single-Progenitor (SP)
% model with balanced or unbalanced cell fates, and experimental data on
% the number of total cells / clone fitted.
% User can select parameter conditions of interest, or choose those used in
% the manuscript (best-fits according to maximum-likelihood inference).
% Properties of total clone sizes are plotted.

%% SELECTION OF MODEL CONDITIONS USED FOR SIMULATIONS:
% selectMODE:
%   'MLE'::     Model and parameter values corresponding to MLE
%   'Custom'::  Set by user along with the parameters below in the code

selectMODE = 'MLE'; % Mode of choice from the ones above

%% DEFINITION OF MODEL PARAMETER CONDITIONS:
% General simulation parameters:
ntime = [0 10.^[-0.8:0.02:1.42]]; % desired time points for data collection (weeks post-induction) - should cover experimental time points!
indiv = 100000; % number of independent simulations (number of simulated clones)

% Choice of model parameter values:
if (strcmpi('MLE',selectMODE))
    % Load default (MLE) parameter values:
    ParamWT = SelectModelParamVal('WT');
    ParamMUT = SelectModelParamVal('Mut');

elseif (strcmpi('Custom',selectMODE))
    % Specific, customized parameter values:
    % WT:
    ParamWT.model = 'SP_balance';
    ParamWT.lambda = 2.9; % division rate (/week)
    ParamWT.tlag = 0.5/7; % refractory period between consecutive divisions (minimum cell-cycle period) (weeks)
    ParamWT.GamShape = 8; % 'shape' parameter of the gamma-distributed cell-cycle period
    ParamWT.dens = 0.52; % proportion of proliferating basal cells
    ParamWT.r = 0.08; % probability of symmetric division leading to two progenitors, PP, or two differentiating cells, DD
    ParamWT.gamma = 3.1417; % stratification rate (/week)
    ParamWT.micetime = [1.4286 4 12 24]; % time points (weeks) when the suprabasal-to-basal cell ratio has been measured
    ParamWT.mice_m_t = [0.7935 1.0417 0.8310 0.5574]; % experimental suprabasal-to-basal cell ratio measurements
    % Mut:
    ParamMUT.model = 'SP_unbalance';
    ParamMUT.lambda = 2.9; % division rate (/week)
    ParamMUT.tlag = 0.5/7; % refractory period between consecutive divisions (minimum cell-cycle period) (weeks)
    ParamMUT.GamShape = 8; % 'shape' parameter of the gamma-distributed cell-cycle period
    ParamMUT.dens = 0.64; % proportion of proliferating basal cells
    ParamMUT.r = 0.17; % default probability of symmetric division leading to two progenitors, PP, or two differentiating cells, DD
    ParamMUT.Delta = 0.03; % fate imbalance favoring PP vs. DD outcome
    ParamMUT.gamma = 5.0734; % stratification rate (/week)
    ParamMUT.micetime = [1.4286 4 12 24]; % time points (weeks) when the suprabasal-to-basal cell ratio has been measured
    ParamMUT.mice_m_t = [0.6014 0.5616 0.5735 0.4430]; % experimental suprabasal-to-basal cell ratio measurements
end

%% COMPUTATIONAL SIMULATION OF TOTAL CLONE DYNAMICS UNDER GIVEN PARAMETERS
% Basal and suprabasal cells are simulated according to the single-progenitor
% (SP) model paradigm under the specific conditions given for WT and Mut
% (either with balanced or unbalanced cell fates, as defined) and total
% clone sizes retrieved.
% The variable suprabasal-to-total cell ratio is considered in the
% simulations, to better estimate the actual shedding rate and thus the
% cellular dynamics in the suprabasal compartment.

% Selection of MonteCarlo simulator (depending on the model type):
% WT:
disp('neutral SP fit for WT clones (balanced fates) | time-dependent shedding rate assumed')
[nx_basal_WT,nx_total_WT,ntime_WT] = MonteCarloSimulator_SP_total_varMu(ntime,ParamWT.lambda,ParamWT.dens,ParamWT.r,ParamWT.gamma,ParamWT.micetime,ParamWT.mice_m_t,indiv,ParamWT.tlag,ParamWT.GamShape);
% Mut:
if strcmpi('SP_balance',ParamMUT.model)
    disp('neutral SP fit for Mut clones (balanced fates) | time-dependent shedding rate assumed')
    [nx_basal_Mut,nx_total_Mut,ntime_Mut] = MonteCarloSimulator_SP_total_varMu(ntime,ParamMUT.lambda,ParamMUT.dens,ParamMUT.r,ParamMUT.gamma,ParamMUT.micetime,ParamMUT.mice_m_t,indiv,ParamMUT.tlag,ParamMUT.GamShape);
elseif strcmpi('SP_unbalance',ParamMUT.model)
    disp('non-neutral SP fit for Mut clones (unbalanced fates) | time-dependent shedding rate assumed')
    [nx_basal_Mut,nx_total_Mut,ntime_Mut] = MonteCarloSimulator_SP_total_unbalance_varMu(ntime,ParamMUT.lambda,ParamMUT.dens,ParamMUT.r,ParamMUT.Delta,ParamMUT.gamma,ParamMUT.micetime,ParamMUT.mice_m_t,indiv,ParamMUT.tlag,ParamMUT.GamShape);
end

% Save data (workspace variables) into file 'myTest.mat' into ./Datasets folder:
if ~exist('Datasets', 'dir') % makes new directory if this does not exist in pwd
   mkdir('Datasets')
end
save ./Datasets/myTest_total.mat

%% CALCULATE AND PLOT AVERAGE TOTAL CLONE SIZE OVER TIME:

% Select whether to calculate and plot plausible intervals on simulation estimates:
showCI = 1; % Show plausible intervals on model outcome ( 0=NO | 1=YES )  (be aware of the significant time consumption)
% (bear in mind plausible intervals reflect the level of uncertainty given a limited experimental sampling, similar to the avg. number of clones per
% experimental time point, reason why they differ from the 95% confidence intervals shown for certain plots in the manuscript, where bounds were given by the
% 95% CI on the model parameters (MLE) and required launching simulations at all those different parameter values - thus not covered in this code for simplicity)

% Plausible interval settings:
if showCI == 1 % A No. of clones aprox. equivalent to the No. of experimental surviving clones at the different time points is sampled for plausible interval calculations
    sampling = struct('NSubsets',20,'NClones',500); % parameters of random permutation sampling for plausible interval calculation
else
    sampling = struct(); % void parameters (plausible interval not requested)
end

% Calculate:
[avgCloneSize_WT.avg,avgCloneSize_WT.ci95up,avgCloneSize_WT.ci95dn] = calculate_AvgCloneSize(ntime_WT,nx_total_WT,nx_basal_WT,showCI,sampling);
[avgCloneSize_Mut.avg,avgCloneSize_Mut.ci95up,avgCloneSize_Mut.ci95dn] = calculate_AvgCloneSize(ntime_Mut,nx_total_Mut,nx_basal_Mut,showCI,sampling);

% Plot (with or without plausible intervals):
FigProp_WT = struct('Color',[0.8 0.8 0],'Leyend','WT','data','Avg. total cells / clone');
plot_AvgCloneSize_simulated(ntime_WT,avgCloneSize_WT,showCI,FigProp_WT)
FigProp_Mut = struct('Color',[0.23 0.44 0.34],'Leyend','Mut','data','Avg. total cells / clone');
plot_AvgCloneSize_simulated(ntime_Mut,avgCloneSize_Mut,showCI,FigProp_Mut)

% Plot experimental data overlaid:
load ./Datasets/LineageTracing_WT_RYFP_dataset.mat
FigProp_WT = struct('Color',[0.8 0.8 0],'YLim',[0 60],'gap',0);
plot_AvgCloneSize_experimental(rtime,nmice,rx_total,rx_total_all,rx_basal_complete,rx_basal_complete_all,1,FigProp_WT)
load ./Datasets/LineageTracing_Mut_RYFP_dataset.mat
FigProp_Mut = struct('Color',[0.23 0.44 0.34],'YLim',[0 60],'gap',1/7);
plot_AvgCloneSize_experimental(rtime,nmice,rx_total,rx_total_all,rx_basal_complete,rx_basal_complete_all,1,FigProp_Mut)

%% CALCULATE AND PLOT TOTAL CLONE SIZE DISTRIBUTIONS OVER TIME:

% Settings:
cutoffB = 2; % minimum number of basal cells for clone size frequency calculations
NSubsets = 100; % number of subsets simulations are partitioned into to build uncertainty bounds on frequency estimates
for aja = 1:length(rtime)
    evalTime(1,aja) = find(ntime_WT>=rtime(aja),1); % selection of simulated time points coincident with experimental ones
end

% Calculate clone size frequencies: WT
load ./Datasets/LineageTracing_WT_RYFP_dataset.mat
% experimental
[rfreq_rel_WT,rfreq_bin_rel_WT,rbin_label_WT,sem_rfreq_rel_WT,sem_rfreq_bin_rel_WT] = calculate_CloneSizeDist(rtime,rx_total_all_raw,rx_basal_complete_all_raw,cutoffB,struct('Type','Experimental','nmice',nmice),[],rx_total,rx_basal_complete);
% simulated
[nfreq_rel_WT,nfreq_bin_rel_WT,nbin_label_WT,nfreq_rel_95ci_WT,nfreq_bin_rel_95ci_WT] = calculate_CloneSizeDist(ntime_WT(evalTime),nx_total_WT(:,evalTime),nx_basal_WT(:,evalTime),cutoffB,struct('Type','Simulated','NSubsets',NSubsets),rx_basal_complete_all_raw);

% Calculate clone size frequencies: Mut
load ./Datasets/LineageTracing_Mut_RYFP_dataset.mat
% experimental
[rfreq_rel_Mut,rfreq_bin_rel_Mut,rbin_label_Mut,sem_rfreq_rel_Mut,sem_rfreq_bin_rel_Mut] = calculate_CloneSizeDist(rtime,rx_total_all_raw,rx_basal_complete_all_raw,cutoffB,struct('Type','Experimental','nmice',nmice),[],rx_total,rx_basal_complete);
% simulated
[nfreq_rel_Mut,nfreq_bin_rel_Mut,nbin_label_Mut,nfreq_rel_95ci_Mut,nfreq_bin_rel_95ci_Mut] = calculate_CloneSizeDist(ntime_Mut(evalTime),nx_total_Mut(:,evalTime),nx_basal_Mut(:,evalTime),cutoffB,struct('Type','Simulated','NSubsets',NSubsets),rx_basal_complete_all_raw);

% Plot distributions of clone size frequencies (with uncertainty bounds):
FigProp_WT = struct('ColorData',[0.7 0.7 0.7],'ColorFit',[0.8 0.8 0],'row',0,'Leyend','WT data','data','Total clone size');
plot_CloneSizeDist(rtime,rfreq_bin_rel_WT,sem_rfreq_bin_rel_WT,nfreq_bin_rel_WT,nfreq_bin_rel_95ci_WT,rbin_label_WT,nbin_label_WT,FigProp_WT);
FigProp_Mut = struct('ColorData',[0.5 0.5 0.5],'ColorFit',[0.23 0.44 0.34],'row',1,'Leyend','Mut data','data','Total clone size');
plot_CloneSizeDist(rtime,rfreq_bin_rel_Mut,sem_rfreq_bin_rel_Mut,nfreq_bin_rel_Mut,nfreq_bin_rel_95ci_Mut,rbin_label_Mut,nbin_label_Mut,FigProp_Mut);


%% COMPUTATIONAL SIMULATION OF TOTAL CLONE DYNAMICS UNDER GIVEN PARAMETERS (ALTERNATIVE)
% In this alternative with time-dependent suprabasal-to-total cell ratio, an
% even induction efficiency over the different cell types is assumed as initial condition.
% The fraction of floating clones at short time scales can be sensitive to the
% initial induction condition, whether this was restricted to proliferating
% basal cells or any cell could be initially labelled. Thus the convenience to
% run simulations where the latter scenario is considered...

% Selection of MonteCarlo simulator (depending on the model type):
% WT:
disp('neutral SP fit for WT clones (balanced fates) | time-dependent shedding rate assumed & even initial induction')
[nx_basal_WT2,nx_total_WT2,ntime_WT2] = MonteCarloSimulator_SP_total_varMu_evenStart(ntime,ParamWT.lambda,ParamWT.dens,ParamWT.r,ParamWT.gamma,ParamWT.micetime,ParamWT.mice_m_t,indiv,ParamWT.tlag,ParamWT.GamShape);
% Mut:
if strcmpi('SP_balance',ParamMUT.model)
    disp('neutral SP fit for Mut clones (balanced fates) | time-dependent shedding rate assumed & even initial induction')
    [nx_basal_Mut2,nx_total_Mut2,ntime_Mut2] = MonteCarloSimulator_SP_total_varMu_evenStart(ntime,ParamMUT.lambda,ParamMUT.dens,ParamMUT.r,ParamMUT.gamma,ParamMUT.micetime,ParamMUT.mice_m_t,indiv,ParamMUT.tlag,ParamMUT.GamShape);
elseif strcmpi('SP_unbalance',ParamMUT.model)
    disp('non-neutral SP fit for Mut clones (unbalanced fates) | time-dependent shedding rate assumed & even initial induction')
    [nx_basal_Mut2,nx_total_Mut2,ntime_Mut2] = MonteCarloSimulator_SP_total_unbalance_varMu_evenStart(ntime,ParamMUT.lambda,ParamMUT.dens,ParamMUT.r,ParamMUT.Delta,ParamMUT.gamma,ParamMUT.micetime,ParamMUT.mice_m_t,indiv,ParamMUT.tlag,ParamMUT.GamShape);
end

%% CALCULATE AND PLOT % FLOATING CLONES (CLONES MISSING BASAL ATTACHMENT) OVER TIME:

% Calculate simulated values: WT
[mean_nfreqFloat_WT] = calculate_FloatingClones(ntime_WT,nx_basal_WT,nx_total_WT,[],struct('Type','Simulated'));
[mean_nfreqFloat_WT2] = calculate_FloatingClones(ntime_WT2,nx_basal_WT2,nx_total_WT2,[],struct('Type','Simulated'));
% Calculate simulated values: Mut
[mean_nfreqFloat_Mut] = calculate_FloatingClones(ntime_Mut,nx_basal_Mut,nx_total_Mut,[],struct('Type','Simulated'));
[mean_nfreqFloat_Mut2] = calculate_FloatingClones(ntime_Mut2,nx_basal_Mut2,nx_total_Mut2,[],struct('Type','Simulated'));

% Plot simulated values: WT
FigProp_WT = struct('Color',[0.8 0.8 0],'Line','-','Leyend','WT');
plot_FloatingClones_simulated(ntime_WT,mean_nfreqFloat_WT,FigProp_WT);
FigProp_WT = struct('Color',[0.8 0.8 0],'Line','--','Leyend','WT (even induc.)');
plot_FloatingClones_simulated(ntime_WT2,mean_nfreqFloat_WT2,FigProp_WT);
% Plot simulated values: Mut
FigProp_Mut = struct('Color',[0.23 0.44 0.34],'Line','-','Leyend','Mut');
plot_FloatingClones_simulated(ntime_Mut,mean_nfreqFloat_Mut,FigProp_Mut);
FigProp_Mut = struct('Color',[0.23 0.44 0.34],'Line','--','Leyend','Mut (even induc.)');
plot_FloatingClones_simulated(ntime_Mut2,mean_nfreqFloat_Mut2,FigProp_Mut);

% Calculate and Plot experimental values overlaid:
load ./Datasets/LineageTracing_WT_RYFP_dataset.mat
[mean_rfreqFloat_WT,sem_rfreqFloat_WT,rfreqFloat_mice_WT] = calculate_FloatingClones(rtime,rx_basal_complete,rx_total,nmice,struct('Type','Experimental'));
FigProp_WT = struct('Color',[0.8 0.8 0],'YLim',[0 60],'gap',0);
plot_FloatingClones_experimental(rtime,mean_rfreqFloat_WT,sem_rfreqFloat_WT,nmice,rfreqFloat_mice_WT,FigProp_WT);
load ./Datasets/LineageTracing_Mut_RYFP_dataset.mat
[mean_rfreqFloat_Mut,sem_rfreqFloat_Mut,rfreqFloat_mice_Mut] = calculate_FloatingClones(rtime,rx_basal_complete,rx_total,nmice,struct('Type','Experimental'));
FigProp_Mut = struct('Color',[0.23 0.44 0.34],'YLim',[0 60],'gap',1/7);
plot_FloatingClones_experimental(rtime,mean_rfreqFloat_Mut,sem_rfreqFloat_Mut,nmice,rfreqFloat_mice_Mut,FigProp_Mut);

