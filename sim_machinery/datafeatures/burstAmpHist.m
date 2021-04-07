function [nb,Apdf] = burstAmpHist(dataX,fsamp,bins,minbs)

XH = abs(hilbert(dataX));
burstinds = SplitVec(find(XH>prctile(XH,75)),'consecutive');
segL = 1000*(cellfun('length',burstinds)/fsamp);

burstinds(segL<minbs) = [];
amp = 1000*(cellfun(@(a) max(XH(a)),burstinds)/fsamp);
[Apdf,nb] =  ksdensity(amp,bins);

% end