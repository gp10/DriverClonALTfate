function [clonesizes_sampled,clonesizes_ref_sampled] = subsampling_simulated_clones(nsamples,sizesamples,cloneSizes,cloneSizes_ref,filter,timepoints)
%% Repeated sampling (random permutation) of simulated clones into subsets: (alternative version)
% It generates subsets with a limited number of randomly-chosen clones from
% the pool of simulated clones by random permutation (with replacement).
% Parameters nsamples and sizesamples specify the total number of subsets
% and the number of clones contained in each, respectively.

% from Herms et al, 2020

%% Input:
% nsamples: number of subsets the data is distributed in.
% sizesamples: number of random clones assigned to each subset.
% clonesizes: matrix of format (:,timepoints) with clone sizes at given time points, used as input
% clonesizes_ref: matrix (:,timepoints) with basal clone sizes at given time points, used to establish the cutoff - same as 'clonesizes' when this refers to basal and not total sizes
% filter: threshold for clone sampling (only clones in which 'clonesizes_ref' >= filter are sampled)
% timepoints: time points when clone sizes were calculated (expressed in weeks)

%% Output:
% clonesizes_sampled: cell array {subsets,timepoints}(:,1) of clone sizes in each subset at each specified time point
% clonesizes_ref_sampled: cell array {subsets,timepoints}(:,1) of basal clone sizes in each subset at each specified time point

%% Sampling (random permutation) of simulated clones:
clonesizes_sampled = {};
clonesizes_ref_sampled = {};

for luptime = 1:size(timepoints,2)
    % Restrict sampling to clones with a number of cells >= filter
    loc_prolif = find(cloneSizes_ref(:,luptime)>=filter);
    rnd_pickpos = [];
    % Subset making:
    for rnd_loop = 1:nsamples
        % random clone sampling:
        rnd_pickpos = loc_prolif(randperm(size(loc_prolif,1),sizesamples(luptime)),1);
        clonesizes_sampled{rnd_loop,luptime} = cloneSizes(rnd_pickpos,luptime);
        clonesizes_ref_sampled{rnd_loop,luptime} = cloneSizes_ref(rnd_pickpos,luptime);
    end
end
