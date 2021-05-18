function [nb,kpdf,lnpdf,segL] = burstDurHist(dataX,fsamp,bins,minbs)

XH = abs(hilbert(dataX));
burstinds = SplitVec(find(XH>prctile(XH,75)),'consecutive');
segL = 1000*(cellfun('length',burstinds)/fsamp);
segL(segL<minbs) = [];
if numel(segL)>2
    [kpdf,nb] =  ksdensity(segL,bins);
else
    kpdf = nan;
    nb = nan;
end

if nargout>2
    if numel(segL)>2
        [raypdf] = fitdist(segL','rayleigh');
        %         [lnpdf] = fitdist(segL','Lognormal');
        [lnpdf] = fitdist(segL','inversegaussian');
    else
        raypdf = nan;
        lnpdf = nan;
        warning('Not enough data to fit a distribution!')
    end
end