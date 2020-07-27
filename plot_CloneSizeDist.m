function plot_CloneSizeDist(timepoints,rfreq_bin_rel,sem_rfreq_bin_rel,nfreq_bin_rel,nfreq_bin_95ci,rbin_label,nbin_label,FigProp)
%% Plots the clone size distributions over time
% Histograms of the number of basal or total cells per clone are displayed
% for the different experimental time points. Experimental frequencies are
% shown with s.e.m. and model fits with 95% plausible intervals built based
% on the actual number of clones counted.

% from Herms et al, 2020

%% Input:
% timepoints: vector of time points (expressed in weeks)
% rfreq_bin_rel: cell array {1,timepoints} of relative frequencies of experimental clone sizes binned in groups of powers of 2
% sem_rfreq_bin_rel: cell array {1,timepoints} with error bounds (s.e.m.) on 'rfreq_bin_rel'
% nfreq_bin_rel: matrix [:,timepoints] of relative frequencies of simulated clone sizes binned in groups of powers of 2
% nfreq_bin_95ci: cell array {1,timepoints} with uncertainty bounds (95% CI) on nfreq_bin_rel
% rbin_label: labels spanning all experimental clone size categories
% nbin_label: labels spanning all simulated clone size categories
% FigProp: structure containing output display properties
    % struct{ColorData==[0.7 0.7 0.7], ColorFit==[0.8 0.8 0], row==0, Leyend=='WT data', data=='Basal clone size'}
        % ColorData: color used for plotting experimental data (bars)
        % ColorFit: color used for plotting simulation fits (lines)
        % row: row-displacement for subplots (0 by default; set higher number to start plotting on a row different to the first one)
        % Leyend: legend for the experimental data condition
        % data: x-axis label

%% Figure panel selection:
% Make handle to a new figure if this does not exist yet:
if isempty(findobj( 'Type', 'Figure', 'Name', 'CloneSizeDistPlot'))
    % makes new figure:
    figure()
    set(gcf,'Name','CloneSizeDistPlot');
else
    % retrieves existing one:
    figure(findobj( 'Type', 'Figure', 'Name', 'CloneSizeDistPlot'))
end

%% Plot overlaid (experimental and simulated) distributions at different times with errorbars/confidence bounds:
myYlim = [0.6 0.5 0.4 0.4];
myYlimTick = [0 0.3 0.6; 0 0.25 0.5; 0 0.2 0.4; 0 0.2 0.4];
if size(rbin_label,2) > size(nbin_label,2); nbin_label = rbin_label; end

for bisp = 1:length(timepoints)
    subplot(3,4,bisp+4*FigProp.row)
    hold on
    barwitherr( sem_rfreq_bin_rel{1,bisp}(3:size(rfreq_bin_rel{1,bisp},1),1),  rfreq_bin_rel{1,bisp}(3:end,1)   ,'FaceColor',FigProp.ColorData,'EdgeColor','none')
    hold on
    plot([1:size(nfreq_bin_rel,1)-2],   nfreq_bin_rel(3:end,bisp),                              '-','Color',FigProp.ColorFit)
    plot([1:size(nfreq_bin_rel,1)-2],   nfreq_bin_95ci{1,bisp}(3:size(nfreq_bin_rel,1),1),      '-', 'Color',FigProp.ColorFit);
    plot([1:size(nfreq_bin_rel,1)-2],   nfreq_bin_95ci{1,bisp}(3:size(nfreq_bin_rel,1),2),      '-', 'Color',FigProp.ColorFit);
    set(gca,'XTick',[1:size(nbin_label,2)-2])
    set(gca,'XTickLabel',nbin_label(3:end))
    rotateXLabels( gca(), 90 )
    title(sprintf('%.1f weeks',timepoints(1,bisp)))
    yl = ylim(gca);
    yl(1) = 0;
    ylim(gca, yl);
    xlim([0 size(nbin_label,2)-1])
    ylim([0 myYlim(bisp)])
    set(gca,'YTick',myYlimTick(bisp,:))
    if FigProp.row == 1
        set(gca,'YDir','reverse')
    end
    grid on
    set(gca,'GridLineStyle',':')
    if bisp == 1; ylabel('Frequency'); end
end
xlabel(FigProp.data)
legend({FigProp.Leyend,'(s.e.m.)','model fit','(95% PI)'})
