function plotModComp_210320(R,cmap,FI)
if nargin<3
    FI = 1;
end
% addpath('C:\Users\Tim\Documents\MATLAB_ADDONS\violin')
% R.plot.confint = 'none';
%     cmap = brewermap(numel(R.modcompplot.NPDsel),'Spectral');
    % cmap = linspecer(R.modcomp.modN);
%% First get probabilities so you can compute model space epsilon
K = 0;
for modID = R.modcomp.modN
    K = K+1;
    dagname = modID{1};
    load([R.rootn 'outputs\' R.out.tag '\' dagname '\modeProbs_' R.out.tag '_' dagname '.mat'])
    A = varo; %i.e. permMod
    if ~isempty(A)
        r2rep = [A.r2rep{:}];
        r2rep(isnan(r2rep) | isinf(r2rep)) = [];
        r2bank{K} = r2rep;
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
K = 0;

for modID = R.modcompplot.NPDsel
    K = K+1;
    dagname = modID{1};
    load([R.rootn 'outputs\' R.out.tag '\' dagname '\modeProbs_' R.out.tag '_' dagname '.mat'])
    A = varo; %i.e. permMod
    avStruc = averageCell(A.par_rep);
    if ~isempty(A)
        r2rep = [A.r2rep{:}];
        r2rep(isnan(r2rep) | isinf(r2rep)) = [];
        pd = fitdist(r2rep','normal');
        x_values = -1:0.01:1;
        y = pdf(pd,x_values);
%         figure(FI)
%         plot(x_values,y,'LineWidth',2)
        
        r2repc = r2rep; % cut-off (for plotting)
        r2repc(r2repc<prctile(r2repc,15)) = []; % cut-off lower outliers (for plotting)
        r2repSave{K} = (r2repc);
        
        %         figure(1)
        %         histogram(r2rep,-1:0.1:1);
        hold on
        
        KL(K) = sum(A.KL(~isnan(A.KL)));
        DKL(K) = sum(A.DKL);
        pmod(K) =sum(r2rep>R.modcomp.modEvi.epspop) / size(r2rep,2);
        
        h =  figure(FI);
        R.plot.cmap = cmap(K,:);
        flag = 0;
        
        list = find([A.r2rep{:}]>R.modcomp.modEvi.epspop);
        if numel(list)>2
            parcat = [];
            
            for i = list
                parcat(:,i) = spm_vec(A.par_rep{i});
            end
            parMean{K} = spm_unvec(mean(parcat,2),A.par_rep{1});
        else
            parMean{K} = A.par_rep{1};
        end
        if ismember(modID,R.modcompplot.NPDsel)
            p = p+1;
            if K == 1; % only plot empirical data once
                empflag = 1;
            else
                empflag = 0;
            end
            hl(p) = R.plot.outFeatCiFx(R,A,empflag,cmap(p,:)); 
        end
        % hl(modID) = plot(1,1);
        if ~flag
            mni = mni +1;
            longlab{mni} = modID;
        end
    else
        r2repSave{K} = nan(size(r2rep));
    end
    
    shortlab{K} = sprintf('M%.f',K);
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

legend(hl,R.modcompplot.NPDsel)

% hl(end+1) = dl;
% longlab{end+1} = 'Data';
figure(FI+20)
subplot(4,1,1)
catA = [];
for i = 1:numel(r2repSave)
catA = [catA repmat(i,1,numel(r2repSave{i}))];
end
violins = violinplot([r2repSave{:}],catA)

for i = 1:numel(R.modcompplot.NPDsel)
violins(i).ViolinColor = cmap(i,:);
end
hold on
plot([0 numel(R.modcompplot.NPDsel)+1],[R.modcomp.modEvi.epspop R.modcomp.modEvi.epspop],'k--')
xlabel('Model'); ylabel('NMRSE'); grid on;
a = gca; a.XTick = 1:numel(R.modcompplot.NPDsel);
a.XTickLabel = shortlab;

subplot(4,1,2)
% TlnK = 2.*log(max(pmod)./pmod);
%  TlnK = log10(pmod); % smallest (closest to zero) is the best
TlnK = -log10(1-pmod); % largest is the best

for i = 1:numel(R.modcompplot.NPDsel)
    %     b = bar(i,-log10(1-pmod(i))); hold on
    b = bar(i,TlnK(i)); hold on
    b.FaceColor = cmap(i,:);
end
a = gca; a.XTick = 1:numel(R.modcompplot.NPDsel); grid on
a.XTickLabel = shortlab;
xlabel('Model'); ylabel('-log_{10} P(M|D)')
xlim([0.5 numel(R.modcompplot.NPDsel)+0.5])

subplot(4,1,3)
for i = 1:numel(R.modcompplot.NPDsel)
    b = bar(i,DKL(i)); hold on
    b.FaceColor = cmap(i,:);
end
a = gca; a.XTick = 1:numel(R.modcompplot.NPDsel);
a.XTickLabel = shortlab;
grid on
xlabel('Model'); ylabel('KL Divergence')
set(gcf,'Position',[277   109   385   895])
xlim([0.5 numel(R.modcompplot.NPDsel)+0.5])
% subplot(3,1,3)
% bar(DKL)
% xlabel('Model'); ylabel('Joint KL Divergence')

subplot(4,1,4)
TlnK = -log10(1-pmod)- log10((DKL./mean(DKL))); % Like a free energy
for i = 1:numel(R.modcompplot.NPDsel)
    b = bar(i,TlnK(i)); hold on
    b.FaceColor = cmap(i,:);
end
a = gca; a.XTick = 1:numel(R.modcompplot.NPDsel);
a.XTickLabel = shortlab;
grid on
xlabel('Model'); ylabel('ACM')
set(gcf,'Position',[277   109   385   895])
xlim([0.5 numel(R.modcompplot.NPDsel)+0.5])



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