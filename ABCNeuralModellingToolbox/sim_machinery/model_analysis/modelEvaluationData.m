function compStruc = modelEvaluationData(R,modelProbs,modID,DKLls,compStruc)
% Collate draws
r2rep = [modelProbs.r2rep{:}];
% Remove failed draws
r2rep(isnan(r2rep) | isinf(r2rep)) = [];

epspop = R.modcomp.modEvi.epspop(R.tmp.confmat(1,modID));

compStruc.r2repSave{modID} = (r2rep);
compStruc.KL(modID) = sum(modelProbs.KL(~isnan(modelProbs.KL)));
compStruc.DKL(modID) = modelProbs.DKL;
compStruc.pmod(modID) = sum(r2rep>epspop)/ size(r2rep,2);


dlist = find(R.tmp.confmat(1,:)==R.tmp.confmat(1,modID));
F1 = -log10(DKLls(modID)./mean(DKLls(dlist)));
F2 = compStruc.pmod;
F2 = -log10(1-F2(modID)); % F2>>0 is more likely
compStruc.ACS(modID) = F1+F2; 


list = find([modelProbs.r2rep{:}]>epspop);
if numel(list)>2
    compStruc.parMean{modID} = averageCell(modelProbs.par_rep);
else
    compStruc.parMean{modID} = modelProbs.par_rep{1};
end
