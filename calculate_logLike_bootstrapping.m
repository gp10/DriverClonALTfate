function [Lbasal_sampled,Lbasal_bin_sampled,Lbasal_t_sampled,Lbasal_bin_t_sampled] = calculate_logLike_bootstrapping(rtime,rx_basal_all,nx_basal,bootstrap,distProp)
%% CALCULATE log-LIKELIHOOD VALUES FOR MODEL SIMULATION ESTIMATES ON EXPERIMENTAL DATA
% Experimental clone sizes are contrasted with simulated clone sizes
% computed under given parameter conditions; distributions of clone sizes
% are built and a log-likelihood value retrieved based on binned frequencies.

% from Herms et al, 2020

%% Input:
% rtime: vector of time points (expressed in weeks)
% rx_basal_all: cell array {1,rtime}(:,1) of experimental clone sizes
% nx_basal: matrix [:,rtime] of simulated clone sizes
% bootstrap: structure containing simulation bootstrapping-related features
    % struct{NSubsets==20, SizeSubsets==290}
        % NSubsets: number of subsets simulations are partitioned into to get robust estimates of the log-Likelihood value
        % SizeSubsets: number of clones sampled per subset (similar to experiments)
% distProp: structure containing experimental-related features
    % struct{freqs_bin_lim=9}
        % freqs_bin_lim: maximum bin category considered for clone sizes (all clones with a size higher than this number are assigned to this category)

%% Output:
% Lbasal_sampled: cell array {s,1} of log-likelihood values for the different bootstrapping subsamples (s = No. of simulated subsamples)
% Lbasal_bin_sampled: cell array {s,1} of log-likelihood values for the different bootstrapping subsamples (s = No. of simulated subsamples) when sizes frequencies are conveniently binnned in groups increasing in powers of two
% Lbasal_t_sampled: cell array {s,1}(1,rtime) of log-likelihood values for the different bootstrapping subsamples, separated per time point (s = No. of simulated subsamples)
% Lbasal_bin_t_sampled: cell array {s,1} of log-likelihood values for the different bootstrapping subsamples, separated per time point (s = No. of simulated subsamples) when sizes frequencies are conveniently binnned in groups increasing in powers of two

%% CLONE SIZE FREQUENCIES OBSERVED IN EXPERIMENTAL DATA:
% COLLECT FREQUENCIES FOR EACH No. OF CELLS, FROM THE EXPERIMENTAL DATA
                            % and
% COLLECT FREQUENCIES IN BINS (POWERS OF 2), FROM THE EXPERIMENTAL DATA:                           
rx_basal_all_sampled = subsampling_experim_clones(bootstrap.NSubsets,bootstrap.SizeSubsets,rx_basal_all,2,rtime);
for boo = 1:bootstrap.NSubsets
    [rfreq_all_sampled{boo,1}, rfreq_all_rel_sampled{boo,1}] = size2freq(rx_basal_all_sampled{boo,1},rtime,1,rx_basal_all_sampled{boo,1},2);
    [rfreq_dim_all_sampled{boo,1}, rfreq_dim_all_rel_sampled{boo,1}, dim_label_sampled{boo,1}] = size2freqbinned(rfreq_all_sampled{boo,1},rx_basal_all_sampled{boo,1},rtime,1);
    rfreq_dim_all_sampled{boo,1} = freq2squeeze(rfreq_dim_all_sampled{boo,1},distProp.freqs_bin_lim,rtime,1);
    rfreq_dim_all_rel_sampled{boo,1} = freq2squeeze(rfreq_dim_all_rel_sampled{boo,1},distProp.freqs_bin_lim,rtime,1);
end

%% CLONE SIZE DISTRIBUTIONS INFERRED BY MODEL SIMULATION:
% GILLESPIE SIMULATION OF THE PDF over time: (done already)
    % and
% COLLECT RELATIVE FREQUENCIES FOR EACH No. OF CELLS, FROM THE SIMULATIONS OF THE MASTER EQUATION:
    % and
% COLLECT FREQUENCIES IN BINS (POWERS OF 2), FROM THE SIMULATIONS OF THE MASTER EQUATION:
[nfreq, nfreq_rel] = size2freq(nx_basal,rtime,2,nx_basal,2);
[nfreq_dim, nfreq_dim_rel, dim_label] = size2freqbinned(nfreq,nx_basal,rtime,2);

%% CALCULATION OF log-LIKELIHOOD VALUE:
Lbasal_t_sampled = {}; Lbasal_sampled = {};
Lbasal_bin_t_sampled = {}; Lbasal_bin_sampled = {};
for boo = 1:bootstrap.NSubsets
    try
        [Lbasal_t_sampled{boo,1}(1,:), Lbasal_sampled{boo,1}] = logLike_calc(rfreq_all_sampled{boo,1},nfreq_rel,rtime);
        [Lbasal_bin_t_sampled{boo,1}(1,:), Lbasal_bin_sampled{boo,1}] = logLike_calc(rfreq_dim_all_sampled{boo,1},nfreq_dim_rel,rtime);
        
    catch
        disp('Could not compute the LogLikelihood!!! Value assigned NaN')
        Lbasal_t_sampled{boo,1}(1,:) = NaN(1,size(rtime,2)); Lbasal_sampled{boo,1} = NaN;
        Lbasal_bin_t_sampled{boo,1}(1,:) = NaN(1,size(rtime,2)); Lbasal_bin_sampled{boo,1} = NaN;
    end
end
