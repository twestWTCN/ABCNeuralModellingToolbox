close all; clear
%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 3- (IV) Multistart Statistics
%%%%%%%%%%%%%%%%%%%%%%%%
%This should link to your repo folder
repopath = 'C:\Users\timot\Documents\GitHub\ABCNeuralModellingToolbox';
% repopath = 'C:\Users\Tim West\Documents\GitHub\ABCNeuralModellingToolbox';
%
%This should be your projectname
projname = 'ABCValidationPaper';
R = ABCAddPaths(repopath,projname);
% R = ABCsetup_partI_STNGPe(R);
% Sets plotting defaults
ABCGraphicsDefaults

R.out.tag = 'figure3_MultiStart'; % Task tag
R.out.dag = sprintf('NPD_STN_GPe_MultiStart_M%.0f',1); % 'All Cross'
load([R.path.rootn '\outputs\' R.path.projectn '\figure3_MultiStart\MSAsave1.mat'],'parConv','N','parSigConv','pMuAct')
format short g

%Subset of parameters (fix while we run more multistarts)
parind = [1 3 7:10];
dataind = [1:N; N+1:2*N];
p = [];
for data = 1:2
X = parConv(parind,dataind(data,:))';
XAct = pMuAct{data}(parind)';
% Xvar = parSigConv(parind,dataind(data,:))';
c = cvpartition(N,'KFold',N);
for i = 1:N
    inSample  = X(~test(c,i),:);
    outSample = X(test(c,i),:);
%     out(i,:) = HZmvntest(inSample)
     T2 = T2Hot1(inSample,0.05,outSample);
     p(i) = T2.p;
     
     T2 = T2Hot1(inSample,0.05,XAct);
     pActCV(i) = T2.p;

end
    
lkhd(data) = sum(p<0.05)/10;
% Test between samples and True
     T2 = T2Hot1(X,0.05,XAct);
     pAct = T2.p;

end
% idx = test(c,1);

% Test betwee


X = parConv';
% Test if the Samples are from the same distribution H0 = Sampes are from
% same Distribution (Szekely & Rizzo energy test)
DepTest2(X(dataind(1,:),:),X(dataind(2,:),:),'test','energy')
% example gives 5x3x2 -> 10 x 4 [G sampleN; ...]
% data is 10x8x2 -> 90 x 2
X = parConv';
X = [repmat(1,N,1) X(dataind(1,:),:); repmat(2,N,1) X(dataind(2,:),:)];
% T2 = T2Hot2ihe(X,0.05);

% MANOVA
figure
X = parConv';
G = [repmat(1,1,N) repmat(2,1,N)];
gplotmatrix(X,[],G,[],'+xo');
[d,p,stats] = manova1(X,G);

c1 = stats.canon(:,1);
c2 = stats.canon(:,2);
figure()
gscatter(c2,c1,G,[],'oxs')
figure()
manovacluster(stats)

% POST HOC ANOVA

X = parConv';
p = [];
for pr = 1:size(X,2)
Y = [X(dataind(1,:),pr) X(dataind(2,:),pr)];
p(pr,:) = anova2(Y,2,'off');
end


% %%  test of variances
%      inSampleVar = var(inSample);
%      outSamplePostVar = Xvar(test(c,i),:);
%      
%      for pt = 1:size(parind,2)
%          F = inSampleVar(pt)./outSamplePostVar(pt);
%          df1 = size(inSample,1); df2 = 1;
%          xunder = 1./max(0,F);
%          ptvar(pt) = 1- fcdf(xunder,df2,df1);
%      end