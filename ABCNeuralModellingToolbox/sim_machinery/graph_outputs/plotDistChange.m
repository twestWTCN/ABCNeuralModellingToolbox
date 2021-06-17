function [] = plotDistChange(R,psave,pPrecSave,pSkewSave,stdev,ii)

for i = 1:2
    if i == 1
        ls = '--';
        nind = 1;
    else
        ls = '-';
        nind = ii;
    end
%{'.params','.noisecov'}
    M = psave(nind).params; M(M==0) = [];
    P = full(pPrecSave{nind}.params); P(P==0) = [];
    S = full(pSkewSave{nind}.params); %S(S==0) = [];
    
    Ma = M(M>-30); 
    P = P(M>-30).*stdev;
    S = S(M>-30);
    
    cmap = linspecer(length(Ma));
    X = -5:.1:5;
    for Q = 1:length(Ma)
        [p,type,coefs] = pearspdf(X,Ma(Q),P(Q)*stdev,S(Q),3);
        plot(X,p,ls,'color',cmap(Q,:))
        hold on
    end
    
end
xlabel('Connection Strength')
ylim([0 1.5]);
xlim(R.SimAn.pOptBound)

