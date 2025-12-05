function par = postDrawMVN(R,Mfit,pOrg,pIndMap,pSigMap,rep)
disp('Drawing from multivariate normal...')
cflag = 0;
x1 = mvnrnd(Mfit.Mu',Mfit.Sigma',rep)';
% x1 = x1(pIndMap,:);
% Plot the mv normal dist
if R.plot.flag == 1
            set(groot,'CurrentFigure',200); 
    clf
    plotProposalDist(R,Mfit,x1,pIndMap,cflag)
end

% setup pars from base
clear base
base = repmat(spm_vec(pOrg),1,rep);
for i = 1:rep
    base(pIndMap,i) = x1(:,i);
    base(pSigMap,i) = base(pSigMap,1);
    par{i} = spm_unvec(base(:,i),pOrg);
end

