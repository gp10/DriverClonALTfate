function plot_FloatingClones_simulated(timepoints,FracFloat,FigProp)
%% Plots the proportion of simulated floating clones over time
% The proportion of simulated clones with only suprabasal cells is plotted
% over time.

% from Herms et al, 2020

%% Input:
% timepoints: vector of time points (expressed in weeks)
% FracFloat: vector [1,timepoints] with the fraction of simulated floating clones over time
% FigProp: structure containing output display properties
    % struct{Line=='-', Color=='r', Leyend=='WT'}
        % Line: line style used for plotting (convenient when plotting different conditions)
        % Color: color used for plotting
        % Leyend: legend

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

%% Plot time course in the % of simulated floating clones:
hold on
plot(timepoints(1,:).*7,FracFloat*100,'LineStyle',FigProp.Line,'Color',FigProp.Color);

%% Legend - retrieve preexisting legend and/or add new entry:
old_legend=findobj(gcf, 'Type', 'Legend');
if isempty(old_legend) % legend is empty; create one
    legend(FigProp.Leyend)
else
    legend([old_legend.String,FigProp.Leyend])
end