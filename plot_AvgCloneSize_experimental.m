function plot_AvgCloneSize_experimental(timepoints,nmice,clonesizes,clonesizes_all,clonesizes_ref,clonesizes_ref_all,showCI,FigProp)
%% Plots the experimental average clone size over time
% The time course of the experimental average (basal or total) clone size
% per mouse is plotted with or without mean ± s.e.m. notches.

% from Herms et al, 2020

%% Input:
% timepoints: vector of time points (expressed in weeks)
% nmice: vector with the No. of mice per time point
% clonesizes: cell array of format {nmice,timepoints}(:,1) containing (basal or total) clone sizes per animal per time point
% clonesizes_all: cell array of format {1,timepoints}(:,1) containing all (basal or total) clone sizes per time point pooled
% clonesizes_ref: cell array of format {nmice,timepoints}(:,1) containing basal clone sizes per animal per time point - same as 'clonesizes' if those refered to basal and not total sizes
% clonesizes_ref_all: cell array of format {1,timepoints}(:,1) containing all basal clone sizes per time point pooled - same as 'clonesizes_all' if those refered to basal and not total sizes
% showCI: include mean ± s.e.m. errorbars ( 0=NO | 1=YES )
% FigProp: structure containing output display properties
    % struct{Color=='r', gap==0.1, YLim==[0 45]}
        % Color: color used for plotting
        % gap: x-axis shift for visual purposes (convenient to avoid overlapping when plotting data from different conditions)
        % YLim: y-axis limits

%% Example:
% timepoints = [1:3];
% nmice = [2 3 1];
% rx_basal{1,1}=ones(5,1); rx_basal{2,1}=ones(7,1).*2; rx_basal{1,2}=ones(4,1).*6; rx_basal{2,2}=ones(11,1).*5; rx_basal{3,2}=ones(5,1).*7; rx_basal{1,3}=ones(3,1).*12;
% rx_basal_all{1,1}=[ones(5,1); ones(7,1).*2]; rx_basal_all{1,2}=[ones(4,1).*6; ones(11,1).*5; ones(5,1).*7]; rx_basal_all{1,3}=[ones(3,1).*12];
% showCI = 1;
% FigProp.Color = 'r'; FigProp.gap = 0; FigProp.YLim = [0 45];
% plot_AvgCloneSize_experimental(timepoints,nmice,rx_basal,rx_basal_all,rx_basal,rx_basal_all,showCI,FigProp);

%% Figure panel selection:
% Make handle to a new figure if this does not exist yet:
if isempty(findobj( 'Type', 'Figure', 'Name', 'avgCloneSizePlot'))
    % makes new figure:
    figure()
    set(gcf,'Name','avgCloneSizePlot');
else
    % retrieves existing one:
    figure(findobj( 'Type', 'Figure', 'Name', 'avgCloneSizePlot'))
end

%% Calculate average size of surviving clones in the different animals and summary statistics:
[mean_clonesizes,sem_clonesizes,mean_indiv_clonesizes] = calculate_AvgCloneSize_experim(timepoints,nmice,clonesizes,clonesizes_all,clonesizes_ref,clonesizes_ref_all);

%% Plot average clone size over time:
hold on
for ata = 1:length(timepoints)
    plot([repmat(timepoints(ata),nmice(ata),1)+FigProp.gap].*7,mean_indiv_clonesizes{1,ata},'o','Color',FigProp.Color,'MarkerFaceColor',FigProp.Color,'MarkerSize',5)
    if showCI == 1 % include mean ± s.e.m. notches
        h = errorbar([timepoints(ata)+FigProp.gap].*7,mean_clonesizes(1,ata),sem_clonesizes(1,ata),'CapSize',0);
        set(h,'Marker','+','Color',FigProp.Color,'LineWidth',0.25)
    end
end
xlabel('Time (days)');
xlim([0 180]); set(gca,'XTick',[0 30 60 90 120 150 180])
ylim(FigProp.YLim)

end

%% Nested function to calculate average basal clone sizes and retrieve summary statistics:
function [mean_Bsize,sem_Bsize,mean_Bsize_indiv,mean_Bsize_All] = calculate_AvgCloneSize_experim(timepoints,nmice,clonesizes,clonesizes_all,clonesizes_ref,clonesizes_ref_all)
    mean_Bsize_indiv = {};
    mean_Bsize = [];
    sem_Bsize = [];
    for ata = 1:size(timepoints,2)
        % mean basal clone size per animal:
        for ete = 1:nmice(ata)
            row_indiv_persis = find(clonesizes_ref{ete,ata} ~= 0);
            mean_Bsize_indiv{1,ata}(ete,1) = mean(clonesizes{ete,ata}(row_indiv_persis,1));
        end
        % Statistics for the different animals (mean, s.e.m.):
        mean_Bsize(1,ata) = mean(mean_Bsize_indiv{1,ata});
        sem_Bsize(1,ata) = std(mean_Bsize_indiv{1,ata},0,1) ./ sqrt(nmice(ata));
        % Statistics in batch (all data from different animals pooled):
        row_All_persis = find(clonesizes_ref_all{:,ata} ~= 0);
        mean_Bsize_All(ata) = mean(clonesizes_all{1,ata}(row_All_persis,1));
    end
end
