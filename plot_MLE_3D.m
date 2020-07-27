function plot_MLE_3D(LikeRatios, mle_mean, Delta_all, Dens_all, r_all)
%% Plots a 3D heatmap of the likelihood ratio-test outcome in the explored parameter space
% Only values within the 95% CI on the MLE are color-coded.

% from Herms et al, 2020

%% Input:
% LikeRatios: matrix [p,m,n] with Likelihood ratio (LR)-test values (p = No. of all values explored for Delta, m = No. of all values explored for rho, n = No. of all values explored for r)
% mle_mean: parameter values corresponding to the max. likelihood estimate (MLE)
% Delta_all: all values explored for Delta, the progenitor fate imbalance
% Dens_all: all values explored for rho, the proportion of proliferating basal cells
% r_all: all values explored for r, the symmetric division prob.

%% Figure panel selection:
% Make handle to a new figure if this does not exist yet:
if isempty(findobj( 'Type', 'Figure', 'Name', 'MLE_GridSearch3DPlot-1'))
    % makes new figure:
    figure()
    set(gcf,'Name','MLE_GridSearch3DPlot-1');
else
    % retrieves existing one:
    figure(findobj( 'Type', 'Figure', 'Name', 'MLE_GridSearch3DPlot-1'))
end

%% Plot 3D heatmap of LR-test values:
mycolmap = colormap(parula);
mycolmap(1,:) = [1 1 1];
colormap(mycolmap)
[X,Y,Z] = meshgrid(r_all,Dens_all,Delta_all);
Sx = [];
Sy = [];
Sz = Delta_all;
h = slice(X, Y, Z, permute(LikeRatios,[2 3 1]), Sx,Sy,Sz);
set(h,'EdgeColor','none');
set(h,'FaceAlpha',0.6);
colorbar
hold on
h2 = contourslice(X, Y, Z, permute(LikeRatios,[2 3 1]), Sx,Sy,Sz);
set(h2,'edgecolor','k','linewidth',.5);
xlabel('Symmetric division prob, \it{r}')
ylabel('Fraction of prolif. cells, \rho')
zlabel('Division fate bias, \Delta');
xlim([0 0.5]); ylim([0 1]); zlim([0 0.14]);
view([-17 8])
grid on;
box on;
camproj('perspective')
hold off

%% Figure panel selection:
% Make handle to a new figure if this does not exist yet:
if isempty(findobj( 'Type', 'Figure', 'Name', 'MLE_GridSearch3DPlot-2'))
    % makes new figure:
    figure()
    set(gcf,'Name','MLE_GridSearch3DPlot-2');
else
    % retrieves existing one:
    figure(findobj( 'Type', 'Figure', 'Name', 'MLE_GridSearch3DPlot-2'))
end

%% Plot 3D heatmap of LR-test values: version with slices for different values of Delta
for aja = 1:length(Delta_all)
    subplot(3,5,aja)
    mycolmap = colormap(parula);
    mycolmap(1,:) = [1 1 1];
    colormap(mycolmap)
    myData(:,:) = LikeRatios(aja,:,:);
    imagesc(myData)
    grid on;
    box on;
    set(gca,'YDir','normal')
    xlabel('\it{r}')
    ylabel('\rho')
    xlim([0.5 41.5])
    ylim([0.5 81.5])
    set(gca,'YTick',[1 21 41 61 81]); set(gca,'YTickLabel',{'0.2','0.4','0.6','0.8','1.0'});
    set(gca,'XTick',[1,11,21,31,41]); set(gca,'XTickLabel',{'0','0.1','0.2','0.3','0.4'});
    title(sprintf('\\Delta = %.2f',Delta_all(aja)))
end

%% Figure panel selection:
% Make handle to a new figure if this does not exist yet:
if isempty(findobj( 'Type', 'Figure', 'Name', 'MLE_GridSearch3DPlot-3'))
    % makes new figure:
    figure()
    set(gcf,'Name','MLE_GridSearch3DPlot-3');
else
    % retrieves existing one:
    figure(findobj( 'Type', 'Figure', 'Name', 'MLE_GridSearch3DPlot-3'))
end

%% Plot 3D heatmap of LR-test values: alternative slice plot
subplot(3,1,1)
%h2 = slice(X, Y, Z, permute(all_Lbasal_filtered,[1 2 3]), Sx,Sy,Sz);
h2 = slice( permute(X,[3 2 1]), permute(Z,[3 2 1]), permute(Y,[3 2 1]), permute(LikeRatios,[1 3 2]), Sx,Sz,Sy);
set(h2,'LineStyle',':') % afects edges
%set(h2,'LineStyle','none') % afects edges
set(h2,'EdgeColor',[0.7 0.7 0.7])
box off
ylim([0 0.07])
set(gca,'YDir','reverse')
xlim([0 0.5])
view([-82 12])
xlabel('Symmetric division prob, \it{r}')
ylabel('Division fate bias, \Delta')
zlabel('Fraction of prolif. cells, \rho')
