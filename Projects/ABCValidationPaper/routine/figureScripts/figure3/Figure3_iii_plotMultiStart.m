close all
R = ABCAddPaths('Rat_NPD','rat_STN_GPe');
R.projectn = 'Rat_NPD'; % Project Name
R.out.tag = 'STN_GPe_ModComp'; % Task tag
R = simannealsetup_NPD_STN_GPe(R);
R.out.dag = sprintf('NPD_STN_GPe_MultiStart_M%.0f',1); % 'All Cross'

load([R.rootn '\routine\rat_STN_GPe\MultiStartAnalysis\MSAsave4.mat'])
format short g
% Check Real Positive
A = [eigvals eigvals/max(abs(eigvals))];

for multiStart = 1:20
    CMDscaled{multiStart} = Y(Inds(1,multiStart):Inds(2,multiStart),:);
end
convMods = find(R2ms>0.01);
convModSel = convMods;
% convModSel(convMods>11) = [];
cmap = brewermap(10,'*Blues');
cmap(11:20,:) = brewermap(10,'*Reds');
i = 0;
figure
subplot(2,2,4)
for multiStart = convMods
    i = i +1;
    %             sc(i) = scatter3(Y(inds,1),Y(inds,2),Y(inds,3),(20.^szvec)*150,cmap(multiStart,:),'.');
    %             hold on
    %             plot3(Y(inds(2:end),1),Y(inds(2:end),2),Y(inds(2:end),3),'color',cmap(multiStart,:))
    
%     sc(i) = scatter(CMDscaled{multiStart}(1:3:end,1),CMDscaled{multiStart}(1:3:end,3),(15.^R2track{multiStart}(1:3:end))*250,cmap(multiStart,:),'.');
    hold on
        p =     plot(CMDscaled{multiStart}(1:3:end,1),CMDscaled{multiStart}(1:3:end,2),'color',cmap(multiStart,:),'LineWidth',1.25);
%     plotVarWidth(CMDscaled{multiStart}(1:3:end,1),CMDscaled{multiStart}(1:3:end,3),2.5.*(15.^R2track{multiStart}(1:3:end)'),cmap(multiStart,:),5)
end
xlabel('Scaling Dimension 1'); ylabel('Scaling Dimension 2');
grid on
title('Multi-start Posteriors on Projected Coordinates')
box on
% legend(sc)

subplot(2,2,3)
% figure
% b = bar(1:numel(pIndMap),mean(parConv,2));
% b.FaceAlpha = 0;
c = bar(1:numel(pIndMap),parConv);
for i = 1:size(c,2)
    c(i).FaceColor = cmap(i,:);
    c(i).EdgeAlpha = 0;
    if i<11
        ig = i;
    legnames{i} = sprintf('%.0fA',ig);
    else
        ig = i-10;
    legnames{i} = sprintf('%.0fB',ig);
    end
end
hold on
b1 = plot((1:numel(pIndMap))-0.25,mean(parConv(:,1:10),2),' ko');
b1.MarkerFaceColor = cmap(1,:);
b1.MarkerEdgeColor = 'none';
e1 = errorbar((1:numel(pIndMap))-0.25,mean(parConv(:,1:10),2),std(parConv(:,1:10),[],2),'.');
e1.Color = cmap(1,:);
e1.LineWidth = 1.2;

b2 = plot((1:numel(pIndMap))+0.25,mean(parConv(:,11:end),2),' ko');
b2.MarkerFaceColor = cmap(11,:);
b2.MarkerEdgeColor = 'none';
e2 = errorbar((1:numel(pIndMap))+0.25,mean(parConv(:,11:end),2),std(parConv(:,11:end),[],2),'.');
e2.Color = cmap(11,:);
e2.LineWidth = 1.2;
g = gca;
g.XTickLabel = parNames;
g.XTickLabelRotation = -45;
grid on
ylabel('Precision Weighted Posteior Means')
title('Multi-start Posterior Means')

% Plot stars
[h p] = ttest(parConv(:,1:10)',parConv(:,11:20)');

X = ((1:numel(pIndMap))).*h;
X(X==0) = NaN;
sc = scatter(X,repmat(1.25,1,numel(pIndMap)),60,'*');
sc.LineWidth = 1.5;
sc.CData = [0 0 0];

leg = legend(c,legnames,'Orientation','vertical');
leg.Box = 'off';
set(leg,'Position',[0.9131 0.1209 0.0803 0.2472]);
set(gcf,'Position',[259         207        1361         615])

for multiStart = convModSel
    % Precision gain
    subplot(2,2,2)
    szvec = R2track{multiStart}; %
    
    y = R2track{multiStart};
    x = (1:size(y,2))';
    expsat = @(B,x) B(1).*(1-exp(-x./B(2))) + B(3);
    opt.MaxIter = 1e4;
    [beta,r,J,Sigma,mse,errorModelInfo,robustw] = nlinfit(x,y',expsat,[0.2 1.5 -1.5],opt);
%     plot(x',y');
%     hold on
%     plot(x,expsat(beta,x))
    convstat(:,multiStart) = beta;
    R2term(multiStart) = R2track{multiStart}(end);
%     plotVarWidth(1:size(parSig{multiStart},2),mean(log10(parSig{multiStart}'),2),2.5.*(15.^szvec),cmap(multiStart,:),5)
    plot(1:size(parSig{multiStart},2),mean(log10(parSig{multiStart}'),2),'color',cmap(multiStart,:),'LineWidth',2);
    hold on
%     scatter(1:size(parSig{multiStart},2),mean(log10(parSig{multiStart}'),2),(15.^szvec)*100,cmap(multiStart,:),'.');
    % Accuracy gain
    subplot(2,2,1)
%     x = ;
%      t = sign(x).*log(1+abs(x)./10^10);
%     plotVarWidth(1:size(R2track{multiStart},2),mean(R2track{multiStart}',2),2.5.*(15.^szvec),cmap(multiStart,:),5)
    plot(1:size(R2track{multiStart},2),mean(R2track{multiStart}',2),'color',cmap(multiStart,:),'LineWidth',2);
    hold on
%     
%     yyaxis right 
%     plot(1:size(kldTrack{multiStart},2),mean(kldTrack{multiStart}',2),'color',cmap(multiStart,:),'LineWidth',2);
%     hold on
%     hold on
%     szvec = mean(parSig{multiStart}',2); %
%     scatter(1:size(R2track{multiStart},2),mean(R2track{multiStart}',2),(15.^szvec)*100,cmap(multiStart,:),'.');
    
end

subplot(2,2,1)
ylabel('R2');    grid on
xlabel('Iteration');
title('Multi-start Model Convergence')
subplot(2,2,2)
ylabel('Log Precision');    grid on
xlabel('Iteration');
title('Multi-start Parameter Inference')
set(gcf,'Position',[72         -14        1237        1009])