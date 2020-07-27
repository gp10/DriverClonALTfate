function plot_MLE_2D(LikeRatios, mle_mean, Dens_all, YScaleFactor, r_all, XScaleFactor, Chi2_95cut_2deg)
%% Plots a 2D heatmap of the likelihood ratio-test outcome in the explored parameter space
% Only values within the 95% CI on the MLE are color-coded.

% from Herms et al, 2020

%% Input:
% LikeRatios: matrix [m,n] with Likelihood ratio (LR)-test values (m = No. of all values explored for rho, n = No. of all values explored for r)
% mle_mean: parameter values corresponding to the max. likelihood estimate (MLE)
% Dens_all: all values explored for rho, the proportion of proliferating basal cells
% YScaleFactor: factor converting rho change steps into heatmap grid unit
% r_all: all values explored for r, the symmetric division prob.
% XScaleFactor: factor converting r change steps into heatmap grid unit
% Chi2_95cut_2deg: LR-test cutoff for building 95% CI on MLE values - based on Chi^2 dist. (2 degrees of freedom) 

%% Figure panel selection:
% Make handle to a new figure if this does not exist yet:
if isempty(findobj( 'Type', 'Figure', 'Name', 'MLE_GridSearch2DPlot'))
    % makes new figure:
    figure()
    set(gcf,'Name','MLE_GridSearch2DPlot');
else
    % retrieves existing one:
    figure(findobj( 'Type', 'Figure', 'Name', 'MLE_GridSearch2DPlot'))
end

%% Fit inferred avg clone size (optional):
[xloc_ALL,yloc_ALL] = find(LikeRatios>=(-2*Chi2_95cut_2deg));
[LineFit,gof] = fit(yloc_ALL,xloc_ALL,'poly1');
scaleLineFit_p1 = LineFit.p1 * (1/YScaleFactor) / (1/XScaleFactor);
scaleLineFit_p2 = LineFit.p2 * (1/YScaleFactor);
Lambda = 2.9;
inferred_inv_est_tau = 1/scaleLineFit_p1 * Lambda;
% Empirical growth rate of the avg clone size:
% est_tau_inv = 0.5361 (0.2101; 0.8621)
% Deconvolution to infer experimental ratio in the Dens_vs_r plot:
experMean_LineFit_p1 = (Lambda/0.5361) * (1/XScaleFactor) / (1/YScaleFactor);
experUpp_LineFit_p1 = (Lambda/0.2101) * (1/XScaleFactor) / (1/YScaleFactor);
experLow_LineFit_p1 = (Lambda/0.8621) * (1/XScaleFactor) / (1/YScaleFactor);
myfun = fittype('k*x+b','dependent',{'y'},'independent',{'x'},'coefficients',{'k', 'b'});
[myfun_fitting,myfun_fitting_stat] = fit(yloc_ALL,xloc_ALL,myfun,'Startpoint',[experMean_LineFit_p1 0],'Lower',[experMean_LineFit_p1 -Inf],'Upper',[experMean_LineFit_p1 100]);

%% Plot 2D heatmap of LR-test values:
hold on
mycolmap = colormap(parula);
mycolmap(1,:) = [1 1 1];
imagesc(LikeRatios)
colormap(mycolmap)
plot([0:102], experMean_LineFit_p1*[0:102] + myfun_fitting.b ,'LineStyle','--','Color',[0.7 0.7 0.7])
plot([0:102], experUpp_LineFit_p1*[0:102] + myfun_fitting.b ,'LineStyle',':','Color',[0.7 0.7 0.7])
plot([0:102], experLow_LineFit_p1*[0:102] + myfun_fitting.b ,'LineStyle',':','Color',[0.7 0.7 0.7])
plot(find(r_all>=mle_mean.r,1), find(Dens_all>=mle_mean.dens,1),'r*')
text(80,80,'N = 3006') % Number of experimental clones counted
set(gca,'YDir','normal')
ylabel('Fraction of prolif. cells, \rho'); xlabel('Symmetric division prob, \it{r}')
xlim([0.5 101.5]);
ylim([-0.5 100.5]);
set(gca,'XTick',[1 21 41 61 81 101]); set(gca,'XTickLabel',{'0','0.1','0.2','0.3','0.4','0.5'});
set(gca,'YTick',[0,20,40,60,80,100]); set(gca,'YTickLabel',{'0','0.2','0.4','0.6','0.8','1.0'});
box on; grid on;
title('WT parameter estimation')
