function plot_FloatingClones_experimental(timepoints,mean_FreqFloat,sem_FreqFloat,nmice,FreqFloat_mice,FigProp)
%% Plots the proportion of experimental floating clones over time
% The proportion of experimental clones with only suprabasal cells is
% plotted over time with error bounds (s.e.m.).

% from Herms et al, 2020

%% Input:
% timepoints: vector of time points (expressed in weeks)
% mean_FreqFloat: vector [1,timepoints] with the fraction of experimental floating clones over time
% sem_FreqFloat: vector [1,timepoints] with the s.e.m. on the fraction of floating clones over time
% nmice: vector [1,timepoints] with the number of mice per time point
% FreqFloat_mice: cell array {1,timepoints} with the fraction of floating clones per mouse at each time
% FigProp: structure containing output display properties
    % struct{Color=='r', gap==0.1, YLim==[0 45]}
        % Color: color used for plotting
        % gap: x-axis shift for visual purposes (convenient to avoid overlapping when plotting data from different conditions)
        % YLim: y-axis limits

%% Figure panel selection:
% Make handle to a new figure if this does not exist yet:
if isempty(findobj( 'Type', 'Figure', 'Name', 'FloatingClonesPlot'))
    % makes new figure:
    figure()
    set(gcf,'Name','FloatingClonesPlot');
else
    % retrieves existing one:
    figure(findobj( 'Type', 'Figure', 'Name', 'FloatingClonesPlot'))
end

%% Plot time course in the % of experimental floating clones:
hold on
for ata = 1:length(timepoints)
    plot([repmat(timepoints(ata),nmice(ata),1)+FigProp.gap].*7,FreqFloat_mice{1,ata}*100,'o','Color',FigProp.Color,'MarkerFaceColor',FigProp.Color,'MarkerSize',5)
    h = errorbar([timepoints(ata)+FigProp.gap].*7,mean_FreqFloat(1,ata)*100,sem_FreqFloat(1,ata)*100,'CapSize',0);
    set(h,'Marker','+','Color',FigProp.Color,'LineWidth',0.25)
end
ylim(FigProp.YLim)
xlim([0 180])
set(gca,'XTick',[0:30:180])
xlabel('Time (days)')
ylabel('% floating clones')
