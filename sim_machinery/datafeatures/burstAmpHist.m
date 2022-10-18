function [nb,Apdf,normpdf,amp] = burstAmpHist(dataX,fsamp,bins,minbs,winmark,dataC)

XH = abs(hilbert(dataX));
burstinds = SplitVec(find(XH>prctile(XH,75)),'consecutive');
if nargin>4
    if ~isempty(winmark)
        burstinds = edgeCorrect(burstinds,winmark);
    end
else
    winmark = [];
end

segL = 1000*(cellfun('length',burstinds)/fsamp);

burstinds(segL<minbs) = [];
amp = (cellfun(@(a) max(XH(a)),burstinds));


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

if nargin>6
    plotBurstConstruction(fsamp,dataC,dataX,XH,burstinds)
    figure
    H = histogram(amp,bins,'Normalization','pdf'); hold on
    H.FaceColor = [131 131 131]/256;
    plot(nb,Apdf,'LineWidth',1.5,'Color',H.FaceColor.*0.4)
    
    a = gca;
        a.YTickLabel = {};
   a.Color = 'none';
   box off; axis square
    xlabel('Amplitude (Z)')
    ylabel('pdf')
end