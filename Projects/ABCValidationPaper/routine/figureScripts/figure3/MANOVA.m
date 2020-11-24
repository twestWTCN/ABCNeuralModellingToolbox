dataind = [1:10; 11:20];
for data = 1:2
% X = parConv(:,1:10)';
X = parConv(:,dataind(data,:))';
c = cvpartition(10,'KFold',10);
for i = 1:10
    inSample = X(~test(c,i),:);
    outSample = X(test(c,i),:);
    out(i,:) = HZmvntest(inSample)
     T2 = T2Hot1(inSample,0.05,outSample);
     p(i) = T2.p;
end
    
lkhd(data) = sum(p<0.05)/10;
end
idx = test(c,1);


X = parConv';
% Test if the Samples are from the same distribution H0 = Sampes are from
% same Distribution (Szekely & Rizzo energy test)
DepTest2(X(dataind(1,:),:),X(dataind(2,:),:),'test','energy')
% example gives 5x3x2 -> 10 x 4 [G sampleN; ...]
% data is 10x8x2 -> 90 x 2
X = parConv';
X = [repmat(1,10,1) X(dataind(1,:),:); repmat(2,10,1) X(dataind(2,:),:)];
% T2 = T2Hot2ihe(X,0.05);

% MANOVA
X = parConv';
G = [repmat(1,1,10) repmat(2,1,10)];
gplotmatrix(X,[],G,[],'+xo');
[d,p,stats] = manova1(X,G);

c1 = stats.canon(:,1);
c2 = stats.canon(:,2);
figure()
gscatter(c2,c1,G,[],'oxs')
figure()
manovacluster(stats)

