function [mle_mean, mle_max95ci, mle_min95ci, maxLike] = calculate_MLE_3D(LogLike_allparam, Delta_all, Dens_all, r_all, Chi2_cutoff_3d)
%% CALCULATES MAXIMUM LIKELIHOOD ESTIMATE (MLE) PARAMETER VALUES WITH CONFIDENCE BOUNDS (3D GRID SEARCH)
% Loads the matrix of log-likelihood values obtained by the 3D parameter grid
% search on {Delta, rho,r} and computes the maximum likelihood estimate (MLE) and the
% confidence interval on the parameter estimates based on likelihood ratio
% (LR) test and the cutoff defined for 3-degrees of freedom.

% from Herms et al, 2020

%% Input:
% LogLike_allparam: matrix [p,m,n] of log-likelihood values
% Delta_all: vector [1,p] with all possible values for Delta, the progenitor fate imbalance
% Dens_all: vector [1,m] with all possible values for rho, the proportion of proliferating basal cells
% r_all: vector [1,n] with all possible values for r, the symmetric division prob.
% Chi2_cutoff_3d: LR-test cutoff (based on Chi^2 distribution)

%% Output:
% mle_mean: structure containing parameter values corresponding to max. likelihood estimate (MLE)
% mle_max95ci: structure containing max. values of the parameters according to the confidence bounds on the MLE
% mle_min95ci: structure containing min. values of the parameters according to the confidence bounds on the MLE
% maxLike: max. log-likelihood value

%% Pre-settings:
if (nargin < 5)
    Chi2_cutoff_3d = 7.82/2;
end

%% OPTIMAL PARAMETER SET (MLE):
% max. log-likelihood value:
maxLike = max(max(max(LogLike_allparam)));
% parameter values corresponding to max. log-likelihood:
[xloc_opti,yloc_opti,zloc_opti] = ind2sub(size(LogLike_allparam),find(LogLike_allparam == maxLike));
mle_mean.Delta = Delta_all(xloc_opti);
mle_mean.dens = Dens_all(yloc_opti);
mle_mean.r = r_all(zloc_opti);

%% CONFIDENCE INTERVALS (BASED ON LIKELIHOOD-RATIO TEST):
all_xloc_min95ci = []; all_xloc_max95ci = [];
all_yloc_min95ci = []; all_yloc_max95ci = [];
all_dens_min95ci = []; all_dens_max95ci = [];
all_r_min95ci = []; all_r_max95ci = [];
Delta_accept = []; Delta_reject = [];

for buz = 1:length(Delta_all) % Delta
    mySubset(:,:) = LogLike_allparam(buz,:,:);
    [xloc_95ci,yloc_95ci] = find(mySubset >= (maxLike - Chi2_cutoff_3d));
    try
        all_xloc_min95ci(1,buz) = min(xloc_95ci); all_xloc_max95ci(1,buz) = max(xloc_95ci); %dens
        all_yloc_min95ci(1,buz) = min(yloc_95ci); all_yloc_max95ci(1,buz) = max(yloc_95ci); %r
        all_dens_min95ci(1,buz) = Dens_all(all_xloc_min95ci(1,buz)); all_dens_max95ci(1,buz) = Dens_all(all_xloc_max95ci(1,buz));
        all_r_min95ci(1,buz) = r_all(all_yloc_min95ci(1,buz)); all_r_max95ci(1,buz) = r_all(all_yloc_max95ci(1,buz));
        Delta_accept = [Delta_accept buz];
    catch % if there are no values within the given ci...
        all_xloc_min95ci(1,buz) = NaN; all_xloc_max95ci(1,buz) = NaN; %dens
        all_yloc_min95ci(1,buz) = NaN; all_yloc_max95ci(1,buz) = NaN; %r
        all_dens_min95ci(1,buz) = NaN; all_dens_max95ci(1,buz) = NaN;
        all_r_min95ci(1,buz) = NaN; all_r_max95ci(1,buz) = NaN;
        Delta_reject = [Delta_reject buz];
    end
end
% Reconstruct the threshold values:
mle_min95ci.Delta = Delta_all(min(Delta_accept)); mle_max95ci.Delta = Delta_all(max(Delta_accept));
mle_min95ci.dens = min(all_dens_min95ci); mle_max95ci.dens = max(all_dens_max95ci);
mle_min95ci.r = min(all_r_min95ci); mle_max95ci.r = max(all_r_max95ci);

% SUMMARY:
mle_mean
mle_min95ci
mle_max95ci
