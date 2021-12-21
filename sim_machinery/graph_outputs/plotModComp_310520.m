function plotModComp_310520(R,cmap,daglist,modnames)
% addpath('C:\Users\Tim\Documents\MATLAB_ADDONS\violin')
% R.plot.confint = 'none';
if nargin<2
    cmap = brewermap(R.modcomp.modN,'Spectral');
    % cmap = linspecer(R.modcomp.modN);
end
%% First get probabilities so you can compute model space epsilon
for modID = 1:numel(R.modcomp.modN)
    if nargin<3
        R.out.dag = sprintf([R.out.tag '_M%.0f'],modID);
    else
        R.out.dag = daglist{modID}
    end
    load([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' R.out.dag '\modeProbs_' R.out.tag '_'  R.out.dag '.mat'])
    permMod = varo; %i.e. permMod
    if ~isempty(permMod)
        ACCrep{modID} = permMod.ACCrep;
        ACCmean(modID) = median(ACCrep{modID});
    end
end

% This gets the joint space epsilon that you can use to compute exceedence
% probability
prct = 50;
[~,id] =  max(ACCmean);
R.modcomp.modEvi.epspop = prctile(ACCrep{id},prct); % threshold becomes median of model fits

p = 0; % plot counter
mni = 0; % Model counter
for modID = numel(R.modcomp.modN):-1:1
    shortlab{modID} = sprintf('M%.f',R.modcomp.modN(modID)); % Make model label
    % Load in the precompute model iterations
    if nargin<3
        R.out.dag = sprintf([R.out.tag '_M%.0f'],modID);
    else
        R.out.dag = daglist{modID};
    end

     [Rout,~,~,~,permMod] = loadABCData_160620(R);   
    R.data = Rout.data;
    
    if ~isempty(permMod)
        % save the model accuracies
        ACCrep = permMod.ACCrep;

        % now get the exceedence probability
        pmod(modID) =sum(ACCrep>R.modcomp.modEvi.epspop) / size(ACCrep,2);
        % parameter divergences
        KL(modID) = sum(permMod.KL(~isnan(permMod.KL))); % sum across all (marginal distributions)
        DKL(modID) = permMod.DKL; % total joint space KL divergence
        MAP(modID) = permMod.MAP;
        ACC{modID} = permMod.ACCrep;
        MSE{modID} = permMod.r2rep;
        
        %% Plot Data Features with Bayesian confidence intervals
        for n = 1:numel(R.data.datatype)
        h(n,1) = subplot(2,2,n);
        end
        flag = 0;

        if ismember(modID,R.modcompplot.NPDsel)
            p = p +1;
            [hl(p), hp, dl, flag] = PlotFeatureConfInt_gen170620(R,permMod,h, cmap(modID,:));
        end
        % hl(modID) = plot(1,1);
        if ~flag
            mni = mni +1;
        end
    else
        r2repSave{modID} = nan(size(ACCrep));
    end
end

figure(1)
    legend([dl hl],['data',modnames(R.modcompplot.NPDsel)]);
set(gcf,'Position',[488.0000  227.4000  611.4000  534.6000]); 
%% Now Plot Results of Model Comparison
%% Combine marginal posterior distribution from model evidence (normalize)
pmodN = (pmod/0.5);
pModDist = (pmodN)./sum(pmodN); % Normalized probability of exceeding the median of the best model

figure(2)
subplot(2,2,1)
MSE = cellfun(@(x) x(~isnan(x) & ~isinf(x) & x>-5),MSE,'UniformOutput',false); % visualization only

% flatten and label for plotting
X = []; G = [];
for i = 1:numel(MSE)
    X = [X MSE{i}];
    G = [G repmat(i,size(MSE{i}))];
end
V = violinplot(X,G);
for i = unique(G)
    V(i).ViolinColor = cmap(i,:);
    V(i).ViolinAlpha = 0.2;
    V(i).ScatterPlot.MarkerFaceAlpha = 0.1;
end
a =gca;
a.XTickLabel = modnames;
a.XTickLabelRotation = 45;
hold on
plot([0 numel(R.modcomp.modN)+1],[R.modcomp.modEvi.epspop R.modcomp.modEvi.epspop],'k--')
ylabel('R2 Error'); 
grid on; box off; axis square
xlim([0.5 numel(R.modcomp.modN)+0.5]);
ylim([0.25 1])

% Now plot probabilities
subplot(2,2,2)
for i = 1:numel(R.modcomp.modN)
        b = bar(i,pmodN(i)); hold on
    b.FaceColor = cmap(i,:);
end
a = gca; a.XTick = 1:numel(R.modcomp.modN); grid on
a.XTickLabel = modnames;
a.XTickLabelRotation = 45;
 ylabel({'Marginal Probability';'P(M|D)'})
xlim([0.5 numel(R.modcomp.modN)+0.5]);
grid on; box off; axis square

subplot(2,2,3)
dklN =(DKL)/sum(DKL);
for i = 1:numel(R.modcomp.modN)
    b = bar(i,dklN(i)); hold on
    b.FaceColor = cmap(i,:);
end
a = gca; a.XTick = 1:numel(R.modcomp.modN);
a.XTickLabel = modnames;
a.XTickLabelRotation = 45;
grid on; box off; axis square
 ylabel({'Normalized ';'Divergence (DKL)'})
xlim([0.5 numel(R.modcomp.modN)+0.5]);% ylim([0 0.15])
ylim([0 0.5])

subplot(2,2,4)
ACS = pModDist -dklN; % + log10(DKL/sum(DKL));
for i = 1:numel(R.modcomp.modN)
    b = bar(i,ACS(i)); hold on
    b.FaceColor = cmap(i,:);
end
a = gca; a.XTick = 1:numel(R.modcomp.modN);
a.XTickLabel = modnames;
a.XTickLabelRotation = 45;
grid on; box off; axis square
xlabel('Model'); ylabel({'Combined Score ';'P(M|D)-DKL'})
xlim([0.5 numel(R.modcomp.modN)+0.5]);
ylim([-0.3 0.5])
    
set(gcf,'Position',[488.0000  227.4000  611.4000  534.6000]); 
