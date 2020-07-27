function [mle_mean, mle_max95ci, mle_min95ci, maxLike] = calculate_MLE_2D(LogLike_allparam, Dens_all, r_all, Chi2_cutoff_2d)
%% CALCULATES MAXIMUM LIKELIHOOD ESTIMATE (MLE) PARAMETER VALUES WITH CONFIDENCE BOUNDS (2D GRID SEARCH)
% Loads the matrix of log-likelihood values obtained by the 2D parameter grid
% search on {rho,r} and computes the maximum likelihood estimate (MLE) and the
% confidence interval on the parameter estimates based on likelihood ratio
% (LR) test and the cutoff defined for 2-degrees of freedom.

% from Herms et al, 2020

%% Input:
% LogLike_allparam: matrix [m,n] of log-likelihood values
% Dens_all: vector [1,m] with all possible values for rho, the proportion of proliferating basal cells
% r_all: vector [1,n] with all possible values for r, the symmetric division prob.
% Chi2_cutoff_2d: LR-test cutoff (based on Chi^2 distribution)

%% Output:
% mle_mean: structure containing parameter values corresponding to max. likelihood estimate (MLE)
% mle_max95ci: structure containing max. values of the parameters according to the confidence bounds on the MLE
% mle_min95ci: structure containing min. values of the parameters according to the confidence bounds on the MLE
% maxLike: max. log-likelihood value

%% Pre-settings:
if (nargin < 4)
    Chi2_cutoff_2d = 5.99/2;
end

%% COMPUTE OPTIMAL PARAMETER SET (MLE):
% max. log-likelihood value:
maxLike = max(max(LogLike_allparam));
% parameter values corresponding to max. log-likelihood:
[xloc_opti,yloc_opti] = ind2sub(size(LogLike_allparam),find(LogLike_allparam == maxLike));
mle_mean.dens = Dens_all(xloc_opti);
mle_mean.r = r_all(yloc_opti);

%% CONFIDENCE INTERVALS (BASED ON LIKELIHOOD-RATIO TEST):
[xloc_95ci,yloc_95ci] = find(LogLike_allparam >= (maxLike - Chi2_cutoff_2d));
try
    xloc_min95ci = min(xloc_95ci); xloc_max95ci = max(xloc_95ci); %dens
    yloc_min95ci = min(yloc_95ci); yloc_max95ci = max(yloc_95ci); %r
    dens_min95ci = Dens_all(xloc_min95ci); dens_max95ci = Dens_all(xloc_max95ci);
    r_min95ci = r_all(yloc_min95ci); r_max95ci = r_all(yloc_max95ci);
catch % if there are no values within the given ci...
    xloc_min95ci = NaN; xloc_max95ci = NaN; %dens
    yloc_min95ci = NaN; yloc_max95ci = NaN; %r
    dens_min95ci = NaN; dens_max95ci = NaN;
    r_min95ci = NaN; r_max95ci = NaN;
end
% Reconstruct the threshold values:
mle_min95ci.dens = dens_min95ci; mle_max95ci.dens = dens_max95ci;
mle_min95ci.r = r_min95ci; mle_max95ci.r = r_max95ci;

% SUMMARY:
mle_mean
mle_min95ci
mle_max95ci
