function [freq_rel,freq_bin_rel,bin_label,freq_rel_errorb,freq_bin_rel_errorb] = calculate_CloneSizeDist(timepoints,cloneSizes,cloneSizes_ref,cutoff,input,refSample,cloneSizes_perMouse,cloneSizes_ref_perMouse)
%% Calculates the clone size frequencies with error bounds for either experimental data or simulations.
% The frequency of clones with certain number of cells is retrieved and
% confidence bounds obtained in the form of s.e.m. (experimental data) or
% 95% CI limits from bootstrapping subsamples (simulated data).

% from Herms et al, 2020

%% Input:
% timepoints: vector of time points (expressed in weeks)
% cloneSizes: cell array (experiments) or matrix (simulations) of size {1,n}(m,1) or [m,n] containing clone sizes (m = No. of clones, n = No. of time points)
% cloneSizes_ref: cell array (experiments) or matrix (simulations) of size {1,n}(m,1) or [m,n] containing corresponding basal clone sizes (m = No. of clones, n = No. of time points) - same as 'cloneSizes' when those refer to basal and not total cells
% cutoff: minimum number of basal cells considered for frequency calculations
% input: structure containing method-specific features
    % struct{Type=='Experimental', nmice==[2,4,3,2], NSubsets==100}
        % Type: string describing the type of input ('Experimental' or 'Simulated'); this conditions the method used for error bound calculation
        % nmice: number of mice used per time point (in the case of Experimental data only)
        % NSubsets: number of subsets simulations are partitioned into to build uncertainty bounds on frequency estimates (in the case of Simulated data only)
% refSample: experimental clone sizes over time used for reference so as to set the same sample sizes (in the case of Simulated data only)
% cloneSizes_perMouse: cell array of size {p,n}(m,1) containing experimental clone sizes (p= No. of mice, m= No. of clones, n= No. of time points) (in the case of Experimental data only)
% cloneSizes_ref_perMouse: cell array of size {p,n}(m,1) containing experimental basal clone sizes (p= No. of mice, m= No. of clones, n= No. of time points) (in the case of Experimental data only) - same as 'cloneSizes_perMouse' when this refers to basal and not total cells

%% Output:
% freq_rel: cell array {1,timepoints}(:,1) or matrix [:,timepoints] containing relative clone size frequencies. Each row r contains relative number of clones with r-1 basal cells (subject to excluding rule -cutoff- above)
% freq_bin_rel: cell array {1,timepoints} or matrix [:,timepoints] containing relative clone size frequencies binned in increasing powers of two. Each row r contains relative number of clones with a number of basal cells in the range 2^(r-3)+1 - 2^(r-2)
% bin_label: labels the number of cells contained per bin.
% freq_rel_errorb: cell array {1,timepoints}(:,1) containing s.e.m. (experimental data) or 95% confidence bound limits (simulated data) on relative clone size frequencies.
% freq_bin_rel_errorb: cell array {1,timepoints}(:,1) containing s.e.m. (experimental data) or 95% confidence bound limits (simulated data) on relative clone size frequencies binned in increasing powers of two.

%% Calculate clone size distributions:
if strcmpi('Experimental',input.Type)
    %% EXPERIMENTAL CLONE SIZE DISTRIBUTIONS:

    % Experimental clone size distributions:
    [freq, freq_rel] = size2freq(cloneSizes,timepoints,1,cloneSizes_ref,cutoff);
    [freq_bin, freq_bin_rel, bin_label] = size2freqbinned(freq,cloneSizes,timepoints,1);
    
    % Experimental clone size distribution errors between mice (s.e.m.):
    [freq_rel_errorb, freq_bin_rel_errorb] = sem_rfreq_mice(timepoints, input.nmice, cloneSizes_perMouse, cloneSizes_ref_perMouse, cutoff);

elseif strcmpi('Simulated',input.Type)
    %% SIMULATED CLONE SIZE DISTRIBUTIONS
    
    % Simulated clone size distributions:
    [freq, freq_rel] = size2freq(cloneSizes,timepoints,2,cloneSizes_ref,cutoff);
    [freq_bin, freq_bin_rel, bin_label] = size2freqbinned(freq,cloneSizes,timepoints,2);
    
    % Simulated clone size distribution uncertainty around MLE (95% CI from bootstrapping samples):
    % (random sampling of the pool of simulations by random permutation to show 95% plausible bounds on the fitting)
    % No. of simulated clones to sample (similar size as the real data)
    [cloneSizes_sampled,cloneSizes_ref_sampled] = sampling_simClones(input.NSubsets,timepoints,refSample,cloneSizes,cloneSizes_ref,cutoff);
    
    % Retrieve basal or total clone size distributions for each sampled set of simulations:
    [all_nfreq_sampled_rel,all_nfreq_sampled_dim_rel] = get_nfreq_sampled(input.NSubsets,timepoints,cloneSizes_sampled,cloneSizes_ref_sampled,cutoff);

    % Reconstruct the lower/upper bounds of freqs for each basal or total clone size:
    [freq_rel_errorb,freq_bin_rel_errorb] = confint_nfreq_sampled(input.NSubsets,timepoints,all_nfreq_sampled_rel,all_nfreq_sampled_dim_rel);
    
end


%% Nested functions:
% Experimental dispersion (sem):
function [sem_mice_rel, sem_mice_dim_rel] = sem_rfreq_mice(rtime,nmice,rx_clones,rx_clones_ref,cutoff)
    sem_mice_rel = {};
    sem_mice_dim_rel = {};
    for iter_time = 1:size(rtime,2)
        rfreq_mice_rel = zeros(500,1); rfreq_mice_dim_rel = zeros(50,1);
        for iter_mice = 1:nmice(iter_time)
            [rfreq_permice, rfreq_permice_rel] = size2freq(rx_clones{iter_mice,iter_time},rtime(iter_time),2,rx_clones_ref{iter_mice,iter_time},cutoff);
            [rfreq_permice_dim, rfreq_permice_dim_rel, dim_label] = size2freqbinned(rfreq_permice,rx_clones{iter_mice,iter_time},rtime(iter_time),2);
            rfreq_mice_rel(1:size(rfreq_permice_rel,1),iter_mice) = rfreq_permice_rel;
            rfreq_mice_dim_rel(1:size(rfreq_permice_dim_rel,1),iter_mice) = rfreq_permice_dim_rel;
        end
        sem_mice_rel{1,iter_time} = std(rfreq_mice_rel,0,2) ./ sqrt(nmice(iter_time));
        sem_mice_dim_rel{1,iter_time} = std(rfreq_mice_dim_rel,0,2) ./ sqrt(nmice(iter_time));
    end
end


% Sampling of simulated clones:
function [cloneSizes_sampled,cloneSizes_ref_sampled] = sampling_simClones(nsamples,rtime,refSample,cloneSizes,cloneSizes_ref,cutoff)
    % Determine experimental number of clones:
    No_sampled_clones = [];
    for iter_time = 1:size(rtime,2)
        No_sampled_clones = [No_sampled_clones size(find(refSample{1,iter_time}>=cutoff),1)];
    end
    % Sampling (random permutation) of simulated clones:
    [cloneSizes_sampled,cloneSizes_ref_sampled] = subsampling_simulated_clones(nsamples,No_sampled_clones,cloneSizes,cloneSizes_ref,cutoff,rtime);
end


% Get clone size frequencies across subsets (subsamples) of simulations:
function [all_nfreq_sampled_rel,all_nfreq_sampled_bin_rel] = get_nfreq_sampled(nsamples,rtime,nx_clones_sampled,nx_clones_sampled_ref,cutoff)
    % Retrieve the clone size distributions for each sampled set of simulations:
    all_nfreq_sampled_rel = {};
    all_nfreq_sampled_bin_rel = {};
    for rnd_loop = 1:nsamples
        for luptime = 1:size(rtime,2)
            [nfreq_sampled, nfreq_sampled_rel] = size2freq(nx_clones_sampled{rnd_loop,luptime},rtime(luptime),2,nx_clones_sampled_ref{rnd_loop,luptime},cutoff);
            [nfreq_sampled_bin, nfreq_sampled_bin_rel, bin_label_sampled] = size2freqbinned(nfreq_sampled,nx_clones_sampled{rnd_loop,luptime},rtime(luptime),2);
            all_nfreq_sampled_rel{rnd_loop,luptime} = nfreq_sampled_rel;
            all_nfreq_sampled_bin_rel{rnd_loop,luptime} = nfreq_sampled_bin_rel;
        end
    end
end


% Simulations dispersion (95% confidence bounds):
function [nfreq_95ci,nfreq_dim_95ci] = confint_nfreq_sampled(nsamples,rtime,all_nfreq_sampled_rel,all_nfreq_sampled_dim_rel)
    % Reconstruct the lower/upper bounds of freqs for each clone size:
    nfreq_95ci = {};
    nfreq_dim_95ci = {};
    for luptime = 1:size(rtime,2)
        nfreq_errorspan = zeros(1000,nsamples);
        nfreq_dim_errorspan = zeros(50,nsamples);
        for rnd_loop = 1:nsamples
            nfreq_errorspan(1:size(all_nfreq_sampled_rel{rnd_loop,luptime},1),rnd_loop) = all_nfreq_sampled_rel{rnd_loop,luptime};
            nfreq_dim_errorspan(1:size(all_nfreq_sampled_dim_rel{rnd_loop,luptime},1),rnd_loop) = all_nfreq_sampled_dim_rel{rnd_loop,luptime};
        end
        nfreq_95ci{1,luptime} = quantile(nfreq_errorspan,[0.025 0.975],2);
        nfreq_dim_95ci{1,luptime} = quantile(nfreq_dim_errorspan,[0.025 0.975],2);
    end
end

end
