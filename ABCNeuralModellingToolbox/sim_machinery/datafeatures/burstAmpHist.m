function [nb,Apdf,normpdf,amp] = burstAmpHist(dataX,fsamp,bins,minbs,winmark)

XH = abs(hilbert(dataX));
burstinds = SplitVec(find(XH>prctile(XH,75)),'consecutive');
if ~isempty(winmark)
burstinds = edgeCorrect(burstinds,winmark);
end


segL = 1000*(cellfun('length',burstinds)/fsamp);

burstinds(segL<minbs) = [];
amp = 1000*(cellfun(@(a) max(XH(a)),burstinds)/fsamp);
if numel(burstinds)>2
    [Apdf,nb] =  ksdensity(amp,bins);
else
    Apdf = nan;
    nb = nan;
end

if nargout>2
    if numel(amp)>2
        [normpdf] = fitdist(amp','normal');
    else
        normpdf = nan;
        warning('Not enough data to fit a distribution!')
    end
end