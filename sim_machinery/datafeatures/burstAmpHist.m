function [nb,Apdf,normpdf] = burstAmpHist(dataX,fsamp,bins,minbs)

XH = abs(hilbert(dataX));
burstinds = SplitVec(find(XH>prctile(XH,75)),'consecutive');
segL = 1000*(cellfun('length',burstinds)/fsamp);

burstinds(segL<minbs) = [];
amp = 1000*(cellfun(@(a) max(XH(a)),burstinds)/fsamp);
[Apdf,nb] =  ksdensity(amp,bins);

if nargout>2
    if numel(segL)>2
        [normpdf] = fitdist(segL','normal');
    else
        normpdf = nan;
        warning('Not enough data to fit a distribution!')
    end
end