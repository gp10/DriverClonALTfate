function [mean_FreqFloat,sem_FreqFloat,FreqFloat_mice] = calculate_FloatingClones(timepoints,basalSize,totalSize,nmice,input)
%% CALCULATE FRACTION OF FLOATING CLONES:
% The fraction of clones containing only suprabasal cells is retrieved and
% confidence bounds (s.e.m.) build when using experimental data as input.

% from Herms et al, 2020

%% Input:
% timepoints: vector of time points (expressed in weeks)
% basalSize: cell array {p,n}(m,1) (experiments) or matrix [m,n](simulations) containing basal clone sizes (p = No. of mice, m = No. of clones, n = No. of time points)
% totalSize: cell array {p,n}(m,1) (experiments) or matrix [m,n] (simulations) containing total clone sizes (p = No. of mice, m = No. of clones, n = No. of time points)
% nmice: number of mice used per time point (in the case of Experimental data only)
% input: structure containing method-specific features
    % struct{Type=='Experimental', nmice==[2,4,3,2], NSubsets==100}
        % Type: string describing the type of input ('Experimental' or 'Simulated'); this conditions the calculations to be done

%% Output:
% mean_FreqFloat: matrix [1,timepoints] containing mean fraction of floating clones per time point
% sem_FreqFloat: matrix [1,timepoints] containing the s.e.m. per time point (for Experimental data only)
% FreqFloat_mice: cell array {1,n}(p,1) containing fraction of floating clones specified per mouse (for Experimental data only) (p = No. of mice, n = No. of time points)

%% Calculate fraction of floating clones:
sem_FreqFloat = [];
FreqFloat_mice = {};

if strcmpi('Experimental',input.Type)
    %% EXPERIMENTAL % FLOATING CLONES:
    
    % Experimental % floating clones:
    [mean_FreqFloat,sem_FreqFloat,FreqFloat_mice] = meanFreqFloat(timepoints,nmice,basalSize,totalSize);

elseif strcmpi('Simulated',input.Type)
    %% SIMULATED % FLOATING CLONES:
    
    % Simulated % floating clones:
    [mean_FreqFloat] = simulFreqFloat(timepoints,basalSize,totalSize);
    
end


%% Nested functions:
% Experimental avg. and dispersion (s.e.m.)
function [mean_FreqFloat,sem_FreqFloat,FreqFloat_indiv] = meanFreqFloat(timepoints,nmice,rx_basal_complete,rx_total)
    for aja = 1:length(timepoints)
        for eje = 1:nmice(aja)
            FreqFloat_indiv{1,aja}(eje,1) = size(find( rx_total{eje,aja}( find(rx_basal_complete{eje,aja}==0) ) ~=0 ),1) ./ size(find(rx_total{eje,aja}~=0),1);
        end
        % Statistics per animal:
        mean_FreqFloat(1,aja) = mean(FreqFloat_indiv{1,aja});
        sem_FreqFloat(1,aja) = std(FreqFloat_indiv{1,aja},0,1) ./ sqrt(nmice(aja));
    end
end

% Simulation avg.
function [PercFloat] = simulFreqFloat(timepoints,nx_basal,nx_total)
    PercFloat = [];
    for ata = 1:length(timepoints)
        PercFloat(ata) = size(find(nx_total(nx_basal(:,ata)==0,ata)~=0),1) ./ size(find(nx_total(:,ata)~=0),1);
    end
end

end
