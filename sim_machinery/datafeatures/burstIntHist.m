function [nb,Apdf,normpdf,bint] = burstIntHist(dataX,fsamp,bins,minbs)

XH = abs(hilbert(dataX));
burstinds = SplitVec(find(XH>prctile(XH,75)),'consecutive');
segL = 1000*(cellfun('length',burstinds)/fsamp);

burstinds(segL<minbs) = [];
bint = [];
for i = 1:numel(burstinds)-1
    bint(i) = burstinds{i+1}(1)-burstinds{i}(end) - 2;
end

if numel(burstinds)>2
    [Apdf,nb] =  ksdensity(bint,bins);
else
    warning('Not enough data to fit a distribution!')
    Apdf = nan;
    nb = nan;
end

% if nargout>2
%     if numel(segL)>2
%         [normpdf] = fitdist(segL','normal');
%     else
%         normpdf = nan;
%         warning('Not enough data to fit a distribution!')
%     end
% end