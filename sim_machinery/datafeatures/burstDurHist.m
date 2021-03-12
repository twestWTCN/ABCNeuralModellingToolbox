function [nb,nbs,E,raypdf,lnpdf,X] = burstDurHist(dataX,band,fsamp,bins)

X = bandpass(dataX,band,fsamp);
XH = abs(hilbert(X));
burstinds = SplitVec(find(XH>prctile(XH,75)),'consecutive');
segL = 1000*(cellfun('length',burstinds)/fsamp);

     nb = ksdensity(segL,bins);
nbs= nb;
E = nan(size(nb));
% [nb,E] = histcounts(segL,bins,'Normalization','pdf');
% nbs = smooth(nb,4);

if nargout>4
    [raypdf] = fitdist(segL,'rayleigh');
    [lnpdf] = fitdist(segL,'Lognormal');
end