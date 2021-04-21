function [nb,kpdf,lnpdf] = burstDurHist(dataX,fsamp,bins,minbs)

XH = abs(hilbert(dataX));
burstinds = SplitVec(find(XH>prctile(XH,75)),'consecutive');
segL = 1000*(cellfun('length',burstinds)/fsamp);
segL(segL<minbs) = [];
[kpdf,nb] =  ksdensity(segL,bins);

if nargout>2
    if numel(segL)>2
        [raypdf] = fitdist(segL','rayleigh');
        [lnpdf] = fitdist(segL','Lognormal');
    else
        raypdf = nan;
        lnpdf = nan;
        warning('Not enough data to fit a distribution!')
    end
end