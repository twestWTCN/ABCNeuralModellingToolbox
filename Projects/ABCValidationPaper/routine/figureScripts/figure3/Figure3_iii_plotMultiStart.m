close all; clear
%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 3- (II) Multistart Analysis
%%%%%%%%%%%%%%%%%%%%%%%%
%This should link to your repo folder
% repopath = 'C:\Users\timot\Documents\GitHub\ABCNeuralModellingToolbox';
repopath = 'C:\Users\Tim West\Documents\GitHub\ABCNeuralModellingToolbox'
% 
%This should be your projectname
projname = 'ABCValidationPaper';
R = ABCAddPaths(repopath,projname);
R = ABCsetup_partI_STNGPe(R);

R.out.tag = 'figure3_MultiStart'; % Task tag
R.out.dag = sprintf('NPD_STN_GPe_MultiStart_M%.0f',1); % 'All Cross'

load([R.path.rootn '\outputs\' R.path.projectn '\MultiStartAnalysis\MSAsave1.mat'])
format short g

N = 3;
% Perform CMD
T = [parWeighted{:}];
D = pdist(T','euclidean');
[Y,eigvals] = cmdscale(squareform(D));

% Check Real Positive
A = [eigvals eigvals/max(abs(eigvals))];

for multiStart = 1:(N*2)
    CMDscaled{multiStart} = Y(Inds(1,multiStart):Inds(2,multiStart),:);
end
convMods = find(R2ms<0.01);
convModSel = convMods;
% convModSel(convMods>11) = [];
cmap = brewermap(N,'*Blues');
cmap(N+1:N*2,:) = brewermap(N,'*Reds');
i = 0;


figure
subplot(2,2,4)
for multiStart = convMods
    CMD_A = CMDscaled{multiStart}(1:3:end,1);
    CMD_B = CMDscaled{multiStart}(1:3:end,3);
    CMD_C = CMDscaled{multiStart}(1:3:end,2);
    cscl = linspace(10,300,size(CMD_A,1));
    i = i +1;
    sc(i) = scatter(CMD_A,CMD_B,cscl,cmap(multiStart,:),'.');
    hold on
    p =     plot(CMD_A,CMD_B,'color',cmap(multiStart,:),'LineWidth',1.25);
%     plot3(CMD_A,CMD_B,CMD_C,'color',cmap(multiStart,:))
hold on
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
c = bar(1:numel(pMuMap),parConv);
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
b1 = plot((1:numel(pMuMap))-0.25,mean(parConv(:,1:N),2),' ko');
b1.MarkerFaceColor = cmap(1,:);
b1.MarkerEdgeColor = 'none';
e1 = errorbar((1:numel(pMuMap))-0.25,mean(parConv(:,1:N),2),std(parConv(:,1:N),[],2),'.');
e1.Color = cmap(1,:);
e1.LineWidth = 1.2;

b2 = plot((1:numel(pMuMap))+0.25,mean(parConv(:,11:end),2),' ko');
b2.MarkerFaceColor = cmap(N+1,:);
b2.MarkerEdgeColor = 'none';
e2 = errorbar((1:numel(pMuMap))+0.25,mean(parConv(:,N+1:2*N),2),std(parConv(:,N+1:2*N),[],2),'.');
e2.Color = cmap(N+1,:);
e2.LineWidth = 1.2;
g = gca;
g.XTickLabel = parNames;
g.XTickLabelRotation = -45;
grid on
ylabel('Precision Weighted Posteior Means')
title('Multi-start Posterior Means')

% Plot stars
[h p] = ttest(parConv(:,1:N)',parConv(:,N+1:2*N)');

X = ((1:numel(pMuMap))).*h;
X(X==0) = NaN;
sc = scatter(X,repmat(1.25,1,numel(pMuMap)),60,'*');
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
%     [beta,r,J,Sigma,mse,errorModelInfo,robustw] = nlinfit(x,y',expsat,[0.2 1.5 -1.5],opt);
    %     plot(x',y');
    %     hold on
    %     plot(x,expsat(beta,x))
%     convstat(:,multiStart) = beta;
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