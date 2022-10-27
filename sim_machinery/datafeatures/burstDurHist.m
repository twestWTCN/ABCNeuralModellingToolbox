function [nb,kpdf,lnpdf,segL] = burstDurHist(dataX,fsamp,bins,minbs,winmark,dataC,epsThresh)
if nargin<7
    epsThresh = 75;
end
XH = abs(hilbert(dataX));
burstinds = SplitVec(find(XH>prctile(XH,epsThresh)),'consecutive');
if nargin>4
    if ~isempty(winmark)
        burstinds = edgeCorrect(burstinds,winmark);
    end
end
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

% % if nargin>6
% %     plotBurstConstruction(fsamp,dataC,dataX,XH,burstinds)
% %     figure
% %     H = histogram(segL,bins,'Normalization','pdf'); hold on
% %     H.FaceColor = [65 138 179]/256;
% %     plot(nb,kpdf,'LineWidth',1.5,'Color',H.FaceColor.*0.4)
% %     
% %     a = gca;
% %     a.YTickLabel = {};
% %     a.Color = 'none';
% %     box off; axis square
% %     xlabel('Duration (ms)')
% %     ylabel('pdf')
% %     xlim([0 800])
% % end