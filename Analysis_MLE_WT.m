%% SCRIPT USED FOR PLOTTING log-LIKELIHOOD VALUES AND DETERMINING MLE ON THE WT MODEL:
% Once a parameter grid search is finished and a collection of
% log-likelihood values has been retrieved extending the analysis in
% "Analysis_BasalCloneDynamics_logLike" to multiple parameter values, in
% this code we plot the log-likelihood values and determine the maximum
% likelihood estimate (MLE) as the one producing the best fit.

%% SETTINGS:
Chi2_95cut_2deg = 5.99/2; % Likelihood Ratio-test cutoff for building 95% CI on MLE values - based on Chi^2 dist. (2 degrees of freedom) 
DoSmoothing = 1; % Perform noise smoothing on log-likelihood values bewteen adjacent parameter values (0=NO | 1=YES)

%% Parameter values considered in the parameter grid search (2D):
Dens_all = 0.01:0.01:1; % possible values for rho, the proportion of proliferating basal cells
r_all = 0:0.005:0.5; % possible values for r, the symmetric division prob.
YScaleFactor = 100; % factor converting rho change steps into heatmap grid unit
XScaleFactor = 200; % factor converting r change steps into heatmap grid unit

%% RETRIEVE MATRIX OF log-LIKELIHOOD VALUES:
% Load previously calculated logLikelihood values (from individual bootstrapping samples):
load ./Datasets/logLikeValues_WT_2DgridSearch_allrho_allr.mat
for buz = 1:size(Lbasal_bin_allrho_sampled,1)
    mySubset = Lbasal_bin_allrho_sampled{buz,1};
    allmySubset(buz,:,:) = mySubset;%
end

% Retrieve mean logLikelihood values from bootstrapping samples:
meanLogLike(:,:) = mean(allmySubset,1);%

%% CALCULATIONS:
% Smoothing neighbour values: (optional)
% (this is for visual purposes, to reduce noise in MLE estimation)
if DoSmoothing
    meanLogLike_s = meanLogLike;
    for aja = 1:length(meanLogLike(:))
        try
            % values at non-edges of parameter grid
            [I,J] = ind2sub([100,101],aja);
            indNeigh = sub2ind([100,101],[I-1 I-1 I-1 I I I I+1 I+1 I+1],[J-1 J J+1 J-1 J J+1 J-1 J J+1]);
            meanLogLike_s(aja) = mean(meanLogLike(indNeigh));
        catch
            % values at edges of parameter grid
            meanLogLike_s(aja) = meanLogLike(aja);
        end
    end
    meanLogLike = meanLogLike_s;
end

% CALCULATE MLE on parameter values (with 95% confidence bounds):
[mle_mean, mle_max95ci, mle_min95ci, maxLike] = calculate_MLE_2D(meanLogLike, Dens_all, r_all, Chi2_95cut_2deg);

% Calculate matrix of Likelihood Ratio (LR)-test statistic:
LikeRatios = (meanLogLike-maxLike).*2;
LikeRatios(find(LikeRatios<(-2*Chi2_95cut_2deg))) = NaN;

%% PLOTTING:
% plot heatmap of Likelihood Ratio-test statistic:
plot_MLE_2D(LikeRatios, mle_mean, Dens_all, YScaleFactor, r_all, XScaleFactor, Chi2_95cut_2deg)

