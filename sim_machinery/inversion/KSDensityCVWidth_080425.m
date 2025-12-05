function u0 = KSDensityCVWidth_080425(x,xf,W,range,n,fx)
% Compute using default
[f1,x1,u1] = ksdensity(x,xf,'function','pdf','Weights',W);
% try other bandwidths
uu = logspace(log10(u1.*10^range(1)),log10(u1.*10^range(2)),n);
uu = max(min(uu, 5 * std(x)), 0.05 * std(x));
v = zeros(size(uu));
% using same partition each time reduces variation 
cp = cvpartition(length(x),'kfold',6);
for j=1:length(uu)
      % compute log likelihood for test data based on training data
      loglik = @(xtr,xte) sum(log(ksdensity(xtr(:,1),xte(:,1),'function','pdf','Weights',xtr(:,2),'width',uu(j))));
      % sum across all train/test partitions
      v(j) = sum(crossval(loglik,[x' W'],'partition',cp));
end 

[~,maxi] = max(v);
u0 = uu(maxi);

% figure(1)
% clf
% subplot(2,1,1)
% [f0,x0,u0] = ksdensity(x,xf,'function','pdf','Weights',W,'width',u0);
% scatter(x1,f1); hold on; scatter(x0,f0)
% histogram(x,25,'Normalization','pdf')
% subplot(2,1,2)
% 
% [fl0,xl0] = ksdensity(x,xf,'function',fx,'Weights',W,'width',u0);
% [fl1,xl1] = ksdensity(x,xf,'function',fx,'Weights',W,'width',u1);
% scatter(xl1,fl1); hold on; scatter(xl0,fl0)
% 

