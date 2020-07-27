function [ParamVal] = SelectModelParamVal(selectCondition)
%% SELECT MODEL PARAMETER VALUES DEPENDING ON CONDITION OF INTEREST:
% A collection of some preset parameter values used in the manuscript to
% simulate different experimental conditions is provided. In particular,
% the MLE parameters providing the best experimental fits are available.
% The function helps selecting the specific set of interest.

% from Herms et al, 2020

%% Input:
% selectCondition: string; type of simulated experimental condition to load parameter values for ('WT' or 'Mut')

%% Output:
% ParamVal: structure containing the parameter values of interest for that simulated condition

%% Parameter selection:
switch selectCondition

    case 'WT'
        %% SUITABLE PARAMETER VALUES FOR: WT (MLE values for a neutral SP model)
        ParamVal.model = 'SP_balance';
        ParamVal.lambda = 2.9; % division rate (/week)
        ParamVal.tlag = 0.5/7; % refractory period between consecutive divisions (minimum cell-cycle period) (weeks)
        ParamVal.GamShape = 8; % 'shape' parameter of the gamma-distributed cell-cycle period
        ParamVal.dens = 0.52; % proportion of proliferating basal cells
        ParamVal.r = 0.08; % probability of symmetric division leading to two progenitors, PP, or two differentiating cells, DD
        ParamVal.gamma = 3.1417; % stratification rate (/week)
        ParamVal.mu = 1.9089; % average shedding rate (/week) (assuming fixed suprabasal-to-basal cell ratio)
        ParamVal.micetime = [1.4286 4 12 24]; % time points (weeks) when the suprabasal-to-basal cell ratio has been measured
        ParamVal.mice_m_t = [0.7935 1.0417 0.8310 0.5574]; % experimental suprabasal-to-basal cell ratio measurements
        
    case 'Mut'
        %% SUITABLE PARAMETER VALUES FOR: Mut (MLE values for a SP model with unbalance)
        ParamVal.model = 'SP_unbalance';
        ParamVal.lambda = 2.9; % division rate (/week)
        ParamVal.tlag = 0.5/7; % refractory period between consecutive divisions (minimum cell-cycle period) (weeks)
        ParamVal.GamShape = 8; % 'shape' parameter of the gamma-distributed cell-cycle period
        ParamVal.dens = 0.64; % proportion of proliferating basal cells
        ParamVal.r = 0.17; % default probability of symmetric division leading to two progenitors, PP, or two differentiating cells, DD
        ParamVal.Delta = 0.03; % fate imbalance favoring PP vs. DD outcome
        ParamVal.gamma = 5.0734; % stratification rate (/week)
        ParamVal.mu = 3.4107; % average shedding rate (/week) (assuming fixed suprabasal-to-basal cell ratio)
        ParamVal.micetime = [1.4286 4 12 24]; % time points (weeks) when the suprabasal-to-basal cell ratio has been measured
        ParamVal.mice_m_t = [0.6014 0.5616 0.5735 0.4430]; % experimental suprabasal-to-basal cell ratio measurements

end
