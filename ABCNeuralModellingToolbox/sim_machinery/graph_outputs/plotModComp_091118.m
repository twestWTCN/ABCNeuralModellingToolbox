function plotModComp_091118(R,cmap)
% addpath('C:\Users\Tim\Documents\MATLAB_ADDONS\violin')
% R.plot.confint = 'none';
if nargin<2
    cmap = brewermap(R.modcomp.modN,'Spectral');
    % cmap = linspecer(R.modcomp.modN);
end
%% First get probabilities so you can compute model space epsilon
for modID = 1:numel(R.modcomp.modN)
    dagname = sprintf([R.out.tag '_M%.0f'],R.modcomp.modN(modID));
    load([R.rootn 'outputs\' R.out.tag '\NPD_' dagname '\modeProbs_' R.out.tag '_NPD_' dagname '.mat'])
    A = varo; %i.e. permMod
    if ~isempty(A)
        r2rep = [A.r2rep{:}];
        r2rep(isnan(r2rep) | isinf(r2rep)) = [];
        r2bank{modID} = r2rep;
    end
end

prct = 50;
% r2med = median(vertcat(r2bank{:}),2);
% find(r2med>prctile(r2med,25))
% r2bankcat = horzcat(r2bank{[1:4 7:10]});
r2bankcat = horzcat(r2bank{:});

R.modcomp.modEvi.epspop = prctile(r2bankcat,prct); % threshold becomes median of model fits
% Adjust the acceptance threshold if any models have no rejections
exc = ones(1,numel(R.modcomp.modN));
while any(exc==1)
    r2bankcat = horzcat(r2bank{:});
    R.modcomp.modEvi.epspop = prctile(r2bankcat,prct); % threshold becomes median of model fits
    for modID = 1:numel(R.modcomp.modN)
        exc(modID) = sum(r2bank{modID}>R.modcomp.modEvi.epspop)/size(r2bank{modID},2);
    end
    prct = prct+1;
end

p = 0;
mni = 0;
for modID = 1:numel(R.modcomp.modN)
    dagname = sprintf([R.out.tag '_M%.0f'],R.modcomp.modN(modID));
    load([R.rootn 'outputs\' R.out.tag '\NPD_' dagname '\modeProbs_' R.out.tag '_NPD_' dagname '.mat'])
    A = varo; %i.e. permMod
    avStruc = averageCell(A.par_rep);
    if ~isempty(A)
        r2rep = [A.r2rep{:}];
        r2rep(isnan(r2rep) | isinf(r2rep)) = [];
        pd = fitdist(r2rep','normal');
        x_values = -1:0.01:1;
        y = pdf(pd,x_values);
        figure(1)
        plot(x_values,y,'LineWidth',2)
        
        r2repc = r2rep; % cut-off (for plotting)
        r2repc(r2repc<prctile(r2repc,15)) = []; % cut-off lower outliers (for plotting)
        r2repSave{modID} = (r2repc);
        
        %         figure(1)
        %         histogram(r2rep,-1:0.1:1);
        hold on
        
        KL(modID) = sum(A.KL(~isnan(A.KL)));
        DKL(modID) = sum(A.DKL);
        pmod(modID) =sum(r2rep>R.modcomp.modEvi.epspop) / size(r2rep,2);
        
        h = figure(10);
        R.plot.cmap = cmap(modID,:);
        flag = 0;
        
        list = find([A.r2rep{:}]>R.modcomp.modEvi.epspop);
        if numel(list)>2
            parcat = [];
            
            for i = list
                parcat(:,i) = spm_vec(A.par_rep{i});
            end
            parMean{modID} = spm_unvec(mean(parcat,2),A.par_rep{1});
        else
            parMean{modID} = A.par_rep{1};
        end
        if ismember(modID,R.modcompplot.NPDsel)
            p = p +1;
            [hl(p), hp, dl, flag] = PlotFeatureConfInt_gen060818(R,A,h);
        end
        % hl(modID) = plot(1,1);
        if ~flag
            mni = mni +1;
            longlab{mni} = sprintf('Model %.f',R.modcomp.modN(modID));
        end
    else
        r2repSave{modID} = nan(size(r2rep));
    end
    
    shortlab{modID} = sprintf('M%.f',R.modcomp.modN(modID));
    %     pause(2)
    %     figure(10)
    %     h = findobj(gca,'Type','line');
    %     legend(h([size(h,1):-2:1 1]),longlab)
    
end
% for p = [2 3 4 7 8 12]
%     subplot(4,4,p); ylim([0 0.7])
% end
% Save the model parameter average
save([R.rootn 'outputs\' R.out.tag '\' R.out.tag '_model_parameter_averages'],'parMean')

set(gcf,'Position',[680   112   976   893])

hl(end+1) = dl;
longlab{end+1} = 'Data';
figure(2)
subplot(4,1,1)
violin(r2repSave,'facecolor',cmap,'medc','k:','xlabel',shortlab,...
    'bw',[0.025 0.1 0.1]); % 0.025 0.2 0.2 0.025 0.025 0.025 0.025 0.2 0.2]);
hold on
plot([0 numel(R.modcomp.modN)+1],[R.modcomp.modEvi.epspop R.modcomp.modEvi.epspop],'k--')
xlabel('Model'); ylabel('NMRSE'); grid on; ylim([-1 0.4])
a = gca; a.XTick = 1:numel(R.modcomp.modN);
a.XTickLabel = shortlab;

figure(2)
subplot(4,1,1)
h = findobj(gca,'Type','line');
% legend(hl,{longlab{[R.modcompplot.NPDsel end]}})

subplot(4,1,2)
% TlnK = 2.*log(max(pmod)./pmod);
%  TlnK = log10(pmod); % smallest (closest to zero) is the best
TlnK = -log10(1-pmod); % largest is the best

for i = 1:numel(R.modcomp.modN)
    %     b = bar(i,-log10(1-pmod(i))); hold on
    b = bar(i,TlnK(i)); hold on
    b.FaceColor = cmap(i,:);
end
a = gca; a.XTick = 1:numel(R.modcomp.modN); grid on
a.XTickLabel = shortlab;
xlabel('Model'); ylabel('-log_{10} P(M|D)')
xlim([0.5 numel(R.modcomp.modN)+0.5])

subplot(4,1,3)
for i = 1:numel(R.modcomp.modN)
    b = bar(i,DKL(i)); hold on
    b.FaceColor = cmap(i,:);
end
a = gca; a.XTick = 1:numel(R.modcomp.modN);
a.XTickLabel = shortlab;
grid on
xlabel('Model'); ylabel('KL Divergence')
set(gcf,'Position',[277   109   385   895])
xlim([0.5 numel(R.modcomp.modN)+0.5])
% subplot(3,1,3)
% bar(DKL)
% xlabel('Model'); ylabel('Joint KL Divergence')

subplot(4,1,4)
TlnK = -log10(1-pmod)- log10((DKL./mean(DKL))); % Like a free energy
for i = 1:numel(R.modcomp.modN)
    b = bar(i,TlnK(i)); hold on
    b.FaceColor = cmap(i,:);
end
a = gca; a.XTick = 1:numel(R.modcomp.modN);
a.XTickLabel = shortlab;
grid on
xlabel('Model'); ylabel('KL Divergence')
set(gcf,'Position',[277   109   385   895])
xlim([0.5 numel(R.modcomp.modN)+0.5])



%% SCRIPT GRAVE
% Adjust the acceptance threshold if any models have no rejections
% exc = ones(1,numel(R.modcomp.modN));
% while any(exc==1)
%     r2bankcat = horzcat(r2bank{:});
%     R.modcomp.modEvi.epspop = prctile(r2bankcat,prct); % threshold becomes median of model fits
%     for modID = 1:numel(R.modcomp.modN)
%         exc(modID) = sum(r2bank{modID}>R.modcomp.modEvi.epspop)/size(r2bank{modID},2);
%     end
%     prct = prct+1;
% end