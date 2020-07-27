function [avgCloneSize,avgCloneSize_ci95up,avgCloneSize_ci95dn] = calculate_AvgCloneSize(timepoints,cloneSizes,cloneSizes_ref,showCI,sampling)
%% Calculates the average surviving clone size over the time series.
% The time course of the average clone size is obtained from the clone sizes
% at given time points.

% from Herms et al, 2020

%% Input:
% timepoints: vector of time points (expressed in weeks)
% cloneSizes: matrix of size [m,n] containing clone sizes (m = No. of clones, n = No. of time points)
% cloneSizes_ref: matrix of size [m,n] containing basal clone sizes (m = No. of clones, n = No. of time points) - same as 'cloneSizes' when those refer to basal and not total cells
% showCI: calculate plausible intervals on model outcome (based on experimental sampling error) ( 0=NO | 1=YES )
% sampling: structure containing sampling properties for CI estimation (ignored if showCI==0)
    % struct{NSubsets==20, NClones==500}
        % NSubsets: number of subsets simulation data are distributed to (random sampling with replacement)
        % NClones: total number of clones contained in each subset

%% Output:
% avgCloneSize: vector of average clone size at different time points
% avgCloneSize_ci95up: upper 95% bound on the average clone size (sampling error)
% avgCloneSize_ci95dn: bottom 95% bound on the average clone size (sampling error)

%% Example:
% timepoints = [1:3];
% cloneSizes(:,1) = [1 1 1 1]'; cloneSizes(:,2) = [5 5 5 5]'; cloneSizes(:,3) = [10 10 10 10]';
% cloneSizes_ref(:,1) = [1 1 1 1]'; cloneSizes_ref(:,2) = [5 5 5 5]'; cloneSizes_ref(:,3) = [10 10 10 10]';
% showCI = 1;
% sampling = struct('NSubsets',20,'NClones',500);
% [avgCloneSize,avgCloneSize_ci95up,avgCloneSize_ci95dn] = calculate_AvgCloneSize(timepoints,cloneSizes,cloneSizes_ref,showCI,sampling);

%% Default parameter values:
% Fill in unset optional values:
if (nargin < 3)
    showCI = 0;
    sampling = struct();
end

avgCloneSize = [];
avgCloneSize_ci95up = [];
avgCloneSize_ci95dn = [];

%% Calculate average size of surviving clones:
for aja = 1:length(timepoints)
    avgCloneSize(aja) = mean(cloneSizes(find(cloneSizes_ref(:,aja)~=0),aja));
end

%% Calculate plausible intervals between simulation runs (by clonal subsampling):
if showCI == 1
    % Subsampling with replacement:
    SizeSubsets = sampling.NClones .* ones(1,length(timepoints)); % e.g. ~500 clones across time points (to ~fit experimental sample size)
    [cloneSizes_2plot_sampled,cloneSizes_ref_2plot_sampled] = subsampling_clones_overtime(sampling.NSubsets,SizeSubsets,cloneSizes,cloneSizes_ref,1,timepoints);
    % Calculate average clone size from subsamples:
    for aja = 1:length(timepoints)
        for eje = 1:sampling.NSubsets
            avgCloneSize_samples(eje,aja) = mean(cloneSizes_2plot_sampled{eje,1}{1,aja}(find(cloneSizes_ref_2plot_sampled{eje,1}{1,aja}~=0)));
        end
    end
    % Calculate 95% plausible intervals on the average clone size:
    avgCloneSize_ci95up = quantile(avgCloneSize_samples,0.975,1);
    avgCloneSize_ci95dn = quantile(avgCloneSize_samples,0.025,1);
end
