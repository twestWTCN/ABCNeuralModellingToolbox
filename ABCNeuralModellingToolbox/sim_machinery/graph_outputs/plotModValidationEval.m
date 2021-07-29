function plotModValidationEval(R,cmap)
% addpath('C:\Users\Tim\Documents\MATLAB_ADDONS\violin')
% R.plot.confint = 'none';
if nargin<2
    cmap = brewermap(R.modcomp.modN,'Spectral');
    % cmap = linspecer(R.modcomp.modN);
end
%% First get probabilities so you can compute model space epsilon
for modID = 1:numel(R.modcomp.modN)
    SimData = R.tmp.confmat(1,modID);
    SimMod = R.tmp.confmat(2,modID);
    dagname = sprintf('NPD_STN_GPe_ConfMat_DataM%.0f_ParM%.0f',SimData,SimMod); % 'All Cross'
    load([R.rootn 'outputs\' R.out.tag '\' dagname '\modeProbs_' R.out.tag '_' dagname '.mat'])
    A = varo; %i.e. permMod
    if ~isempty(A)
        r2rep = [A.r2rep{:}];
        r2rep(isnan(r2rep) | isinf(r2rep)) = [];
        r2bank{modID} = r2rep;
        DKLls(modID) = sum(A.DKL);
    end
end

% Adjust the acceptance threshold if any models have no rejections
exc = R.tmp.confmat(1,:); % Group by data
R = findComparisonEps(R,r2bank,exc);

p = [0; 0];
compStruc = [];
for modID = 1:numel(R.modcomp.modN)
    SimData = R.tmp.confmat(1,modID);
    SimMod = R.tmp.confmat(2,modID);
    dagname = sprintf('NPD_STN_GPe_ConfMat_DataM%.0f_ParM%.0f',SimData,SimMod); % 'All Cross'
    load([R.rootn 'outputs\' R.out.tag '\' dagname '\modeProbs_' R.out.tag '_' dagname '.mat'])
    modelProbs = varo; %i.e. permMod
    % Short Label
    shortlab{modID} = sprintf('M%.f',R.modcomp.modN(modID));
    
    %
    if ~isempty(modelProbs)
        compStruc = modelEvaluationData(R,modelProbs,modID,DKLls,compStruc);
    else
        compStruc.r2repSave{modID} = nan(size(r2rep));
    end
    
    % Plot the data features
%     flag = 0;
%     h = figure(10);
%     if ismember(modID,R.modcompplot.NPDsel)
%         p(1) = p(1) +1;
%         R.plot.cmap = cmap(modID,:);
%         [hl(p(1)), hp, dl, flag] = PlotFeatureConfInt_gen060818(R,modelProbs,h);
%         if ~flag
%             p(2) = p(2) +1;
%             longlab{p(2)} = sprintf('Model %.f',R.modcomp.modN(modID));
%         end
%     end
end

% Save the model parameter average
save([R.rootn 'outputs\' R.out.tag '\' R.out.tag '_modelCompStruc'],'compStruc')

% Return the original mv colormap
R.plot.cmap = cmap;
% Violin Plots of Accuracy
% hl(end+1) = dl;
% longlab{end+1} = 'Data';
figure(2)
subplot(4,1,1)
plotObjectiveViolin(R,compStruc,shortlab)
% h = findobj(gca,'Type','line');
% % legend(hl,{longlab{[R.modcompplot.NPDsel end]}})
subplot(4,1,2)
plotModelComparisonBars(R,compStruc,shortlab,1)
subplot(4,1,3)
plotModelComparisonBars(R,compStruc,shortlab,2)
subplot(4,1,4)
plotModelComparisonBars(R,compStruc,shortlab,3)


figure
R.plot.cmap = brewermap(100,'*RdYlBu');
subplot(1,3,1)
plotModelComparisonMatrix(R,compStruc,shortlab,1)
subplot(1,3,2)
plotModelComparisonMatrix(R,compStruc,shortlab,2)
subplot(1,3,3)
plotModelComparisonMatrix(R,compStruc,shortlab,3)

% Set figure size
set(gcf,'Position',[294         461        1510         377])


%% SCRIPT GRAVE
