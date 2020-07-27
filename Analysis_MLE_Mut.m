%% SCRIPT USED FOR PLOTTING log-LIKELIHOOD VALUES AND DETERMINING MLE ON THE MUTANT MODEL:
% Once a parameter grid search is finished and a collection of
% log-likelihood values has been retrieved extending the analysis in
% "Analysis_BasalCloneDynamics_logLike" to multiple parameter values, in
% this code we plot the log-likelihood values and determine the maximum
% likelihood estimate (MLE) as the one producing the best fit.

%% SETTINGS:
Chi2_95cut_3deg = 7.82/2; % Likelihood Ratio-test cutoff for building 95% CI on MLE values - based on Chi^2 dist. (3 degrees of freedom) 
DoSmoothing = 0; % Perform noise smoothing on log-likelihood values bewteen adjacent parameter values (0=NO | 1=YES)

%% Parameter values considered in the parameter grid search (3D):
Delta_all = 0:0.01:0.14; % possible values for Delta, the progenitor fate imbalance
Dens_all = 0.2:0.01:0.99; % possible values for rho, the proportion of proliferating basal cells
r_all = 0:0.01:0.40; % possible values for r, the symmetric division prob.
XScaleFactor = 100; % factor converting Delta change steps into heatmap grid unit
ZScaleFactor = 100; % factor converting rho change steps into heatmap grid unit
YScaleFactor = 100; % factor converting r change steps into heatmap grid unit

%% RETRIEVE MATRIX OF log-LIKELIHOOD VALUES:
% Load previously calculated logLikelihood values (from individual bootstrapping samples):
load ./Datasets/logLikeValues_Mut_3DgridSearch_allDelta_allrho_allr.mat
for buz = 1:size(Lbasal_bin_allr_sampled,1)
    mySubset = Lbasal_bin_allr_sampled{buz,1};
    allmySubset(buz,:,:,:) = mySubset;
end

% Retrieve mean logLikelihood values from bootstrapping samples:
meanLogLike(:,:,:) = mean(allmySubset,1);%

%% CALCULATIONS:
% Smoothing neighbour values: (optional)
% (this is for visual purposes, to reduce noise in MLE estimation)
if DoSmoothing
    meanLogLike_s = meanLogLike;
    for is = 1:length(Delta_all)
        for aja = 1:length(meanLogLike(is,:))
            try
                [I,J] = ind2sub([size(meanLogLike,2),size(meanLogLike,3)],aja);
                indNeigh = sub2ind([size(meanLogLike,2),size(meanLogLike,3)],[I-1 I-1 I-1 I I I I+1 I+1 I+1],[J-1 J J+1 J-1 J J+1 J-1 J J+1]);
                %mytarget(:) = Lbasal_bin_allparam(is,indNeigh);
                %mytarget2 = mytarget(find(isnan(mytarget)~=1));
                meanLogLike_s(is,aja) = mean(meanLogLike(is,indNeigh));
                %Lbasal_bin_allparam_s(is,aja) = mean(mytarget2);
            catch
                meanLogLike_s(is,aja) = meanLogLike(is,aja);
            end
        end
    end
    meanLogLike = meanLogLike_s;
end

% CALCULATE MLE on parameter values (with 95% confidence bounds):
[mle_mean, mle_max95ci, mle_min95ci, maxLike] = calculate_MLE_3D(meanLogLike, Delta_all, Dens_all, r_all, Chi2_95cut_3deg);
Lambda = 2.9;
mle_mean.gamma = (mle_mean.dens*Lambda - 2*mle_mean.Delta*mle_mean.r*Lambda) / (1-mle_mean.dens);

% Calculate matrix of Likelihood Ratio (LR)-test statistic:
LikeRatios = (meanLogLike-maxLike).*2;
LikeRatios(find(LikeRatios<(-2*Chi2_95cut_3deg))) = NaN;

%% PLOTTING:
% plot heatmap of Likelihood Ratio-test statistic:
plot_MLE_3D(LikeRatios, mle_mean, Delta_all, Dens_all, r_all)

