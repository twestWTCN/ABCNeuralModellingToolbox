function Figure3_iii_plotMultiStart(R)
%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 3- (III) Multistart Plotting
%%%%%%%%%%%%%%%%%%%%%%%%
% Sets plotting defaults
ABCGraphicsDefaults

R.out.tag = 'figure3_MultiStart'; % Task tag
R.out.dag = sprintf('NPD_STN_GPe_MultiStart_M%.0f',1); % 'All Cross'
load([R.path.rootn '\outputs\' R.path.projectn '\figure3_MultiStart\MSAsave1.mat'])
format short g

% Perform CMD
T = [parMS{:} pMuAct{1} pMuAct{2}];
D = pdist(T','euclidean');
[Y,eigvals] = cmdscale(squareform(D));
% Zero centr
Y = Y - Y(1,:);
% Check Real Positive
A = [eigvals eigvals/max(abs(eigvals))];

% MDS for posteriors
for multiStart = 1:(N*2)
    CMDscaled{multiStart} = Y(Inds(1,multiStart):Inds(2,multiStart),:);
end
% MDS for priors
CMDprior{1} = Y(end-1,:);
CMDprior{2} = Y(end,:);

convMods = 1:(N*2);
cmap = brewermap(N,'*Blues');
cmap(N+1:N*2,:) = brewermap(N,'*Reds');
i = 0;

%% Plot multistart traces
figure
subplot(2,2,4)
for multiStart = convMods
    CMD_A = CMDscaled{multiStart}(1:2:end,1);
    CMD_B = CMDscaled{multiStart}(1:2:end,2);
    CMD_C = CMDscaled{multiStart}(1:2:end,3);
    cscl = linspace(10,300,size(CMD_A,1));
    i = i +1;
    sc(i) = scatter(CMD_A,CMD_B,cscl,cmap(multiStart,:),'.');
    hold on
    p =     plot(CMD_A,CMD_B,'color',cmap(multiStart,:),'LineWidth',1.25);
    %     plot3(CMD_A,CMD_B,CMD_C,'color',cmap(multiStart,:))
    hold on
end

% Plot Priors
scatter(CMDprior{1}(1),CMDprior{1}(2),200,'MarkerEdgeColor',cmap(1,:),'Marker','x','LineWidth',1.25);
scatter(CMDprior{2}(1),CMDprior{2}(2),200,'MarkerEdgeColor',cmap(N+1,:),'Marker','x','LineWidth',1.25);

% xlim([-0.15 1.2]); ylim([-1 0])
xlabel('Scaling Dimension 1'); ylabel('Scaling Dimension 2');
grid on
title('Multi-start Posteriors on Projected Coordinates')
box on

%% Bar plots of parameters
subplot(2,2,3)
c = bar(1:numel(pMuMap),parConv);
for i = 1:size(c,2)
    c(i).FaceColor = cmap(i,:);
    c(i).EdgeAlpha = 0;
    if i<9
        ig = i;
        legnames{i} = sprintf('%.0fA',ig);
    else
        ig = i-8;
        legnames{i} = sprintf('%.0fB',ig);
    end
end
hold on
p = scatter((1:numel(pMuMap))-0.25,pMuAct{1},'x','LineWidth',1.25,'MarkerEdgeColor',cmap(1,:));
p = scatter((1:numel(pMuMap))+0.25,pMuAct{2},'x','LineWidth',1.25,'MarkerEdgeColor',cmap(N+1,:));

hold on
b1 = plot((1:numel(pMuMap))-0.25,mean(parConv(:,1:N),2),' ko');
b1.MarkerFaceColor = cmap(1,:);
b1.MarkerEdgeColor = 'none';
e1 = errorbar((1:numel(pMuMap))-0.25,mean(parConv(:,1:N),2),std(parConv(:,1:N),[],2),'.');
e1.Color = cmap(1,:);
e1.LineWidth = 1.2;

b2 = plot((1:numel(pMuMap))+0.25,mean(parConv(:,N+1:end),2),' ko');
b2.MarkerFaceColor = cmap(N+1,:);
b2.MarkerEdgeColor = 'none';
e2 = errorbar((1:numel(pMuMap))+0.25,mean(parConv(:,N+1:2*N),2),std(parConv(:,N+1:2*N),[],2),'.');
e2.Color = cmap(N+1,:);
e2.LineWidth = 1.2;
g = gca;
g.XTickLabel = parNames;
g.XTickLabelRotation = -45;
grid on; box off
ylabel('Multistart MAP Estimates')
title('Log Scaling Factor')

% Plot statistics stars
[h p] = ttest(parConv(:,1:N)',parConv(:,N+1:2*N)');

X = ((1:numel(pMuMap))).*h;
X(X==0) = NaN;
sc = scatter(X,repmat(0.75,1,numel(pMuMap)),60,'*');
sc.LineWidth = 1.5;
sc.CData = [0 0 0];
 xlim([0.25 10.50]); ylim([-1.5 1])
leg = legend(c,legnames,'Orientation','vertical');
leg.Box = 'off';
set(leg,'Position',[0.9131 0.1209 0.0803 0.2472]);

%% Track the accuracy and estimated parameter precision
for multiStart = convMods
    % Precision gain
    subplot(2,2,2)
    szvec = R2track{multiStart}; %
    
    R2term(multiStart) = R2track{multiStart}(end);
    iprec(multiStart) = mean(log10(parSig{multiStart}(1:end-2,end)'));
        iprecPrior(multiStart) = mean(log10(parSig{multiStart}(1:end-2,1)'));

    plot(1:size(parSig{multiStart},2),mean(log10(parSig{multiStart}(1:end-2,:)'),2),'color',cmap(multiStart,:),'LineWidth',2);
    hold on
    % Accuracy gain
    subplot(2,2,1)
    A = R2track{multiStart}';
    A = -log(-A);
%         yyaxis left
    plot(1:size(R2track{multiStart},2),A,'color',cmap(multiStart,:),'Marker','none','LineStyle','-','LineWidth',2);
    hold on
%     yyaxis right
%     plot(1:size(R2track{multiStart},2),mean(pwpe{multiStart},1),'Marker','none','LineStyle','--','color',cmap(multiStart,:),'LineWidth',2);
end

subplot(2,2,1)
ylabel('R2');    grid on; box on
xlabel('Iteration');
title('MDS Projection of Multistart Parameter Trajectories')
subplot(2,2,2)
ylabel('Log mean Precision');    grid on; box on
xlabel('Iteration');
title('Multi-start Parameter Inference')
set(gcf,'Position',[72         -14        1237        1009])
a = 1;