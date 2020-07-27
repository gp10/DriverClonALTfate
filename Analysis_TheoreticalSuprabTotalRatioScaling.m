%% SCRIPT USED FOR THEORETICAL INFERENCE ON THE EFFECT OF UNBALANCED FATES ON THE SB-to-TOTAL CELL RATIO:

%% LOAD REFERENCE NEUTRAL MODEL (WT MLE PARAMETER CONDITIONS; BALANCED FATES):
% This constitutes the starting point for the inference (Delta*r = 0)
ParamRef = SelectModelParamVal('WT');

%% SIMULATE EFFECT OF A FATE IMBALANCE (Delta*r):
% Range of possible values for the imbalance (i.e. product Delta*r, bounded between 0 and 0.5)
Deltar = [0:0.001:0.5];
SBtoTotal = [];
for aja = 1:length(Deltar)
    dens_mut = (2*Deltar(aja)*ParamRef.lambda + ParamRef.gamma) / (ParamRef.lambda + ParamRef.gamma);
    SBtoTotal(aja) = ParamRef.gamma*(1-dens_mut) / (ParamRef.lambda*dens_mut + ParamRef.mu);
end

%% Figure panel selection:
% Make handle to a new figure if this does not exist yet:
if isempty(findobj( 'Type', 'Figure', 'Name', 'SBTotalRatioScalingPlot'))
    % makes new figure:
    figure()
    set(gcf,'Name','SBTotalRatioScalingPlot');
else
    % retrieves existing one:
    figure(findobj( 'Type', 'Figure', 'Name', 'SBTotalRatioScalingPlot'))
end

%% Plot inferred scaling:
plot(Deltar,SBtoTotal,'Color',[0.23 0.44 0.34])
ylim([0 0.5])
xlabel('progenitor fate imbalance (\Delta*r)')
ylabel('1^{st} suprabasal-to-total cell ratio')