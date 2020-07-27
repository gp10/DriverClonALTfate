function [clonesizes_sampled] = subsampling_experim_clones(nsamples,sizesamples,clonesizes,filter,timepoints)
%% Repeated sampling (random permutation) of experimental clones into subsets:
% It generates subsets with a limited number of randomly-chosen clones from
% the pool of experimental clones by random permutation (with replacement).
% Parameters nsamples and sizesamples specify the total number of subsets
% and the number of clones contained in each, respectively.

% from Herms et al, 2020

%% Input:
% nsamples: number of subsets the data is distributed in.
% sizesamples: number of random clones assigned to each subset.
% clonesizes: cell array {1,timepoints} with clone sizes at given time points, used as input
% filter: threshold for clone sampling (only clones in which 'clonesizes_ref' >= filter are sampled)
% timepoints: time points when clone sizes were calculated (expressed in weeks)

%% Output:
% clonesizes_sampled: cell array {subsets,1}{1,timepoints}(:,1) of clone sizes in each subset at each specified time point

%% Sampling (random permutation) of experimental clones:
clonesizes_sampled = {};

for luptime = 1:size(timepoints,2)
    % Restrict sampling to clones with a number of cells >= filter
    loc_prolif = find(clonesizes{1,luptime}(:,1)>=filter);
    rnd_pickpos = [];
    % Subset making:
    for rnd_loop = 1:nsamples
        % random clone sampling:
        mypermut = randperm(size(loc_prolif,1));
        rnd_pickpos = loc_prolif(mypermut(1:sizesamples),1);
        clonesizes_sampled{rnd_loop,1}{1,luptime} = clonesizes{1,luptime}(rnd_pickpos,1);
    end
end
