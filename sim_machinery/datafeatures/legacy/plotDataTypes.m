function plotDataTypes(ftdata)
L = 0;
cmap = hsv(numel(unique(ftdata.hdr.chantype)));
for ctype = unique(ftdata.hdr.chantype)'
    L = L + 1;
    chind = find(strncmp(ftdata.hdr.chantype,ctype{1},3))
     asa{L} = plot(ftdata.time{25},ftdata.trial{25}(chind,:)'- repmat(L*5,1,numel(chind)),'color',cmap(L,:))
     hold on
     
     setL(L) = asa{L}(1);
end

legend(setL, unique(ftdata.hdr.chantype)')

set(gcf,'Position',[ 488.0000  364.2000  537.8000  397.8000])
%
xind = strncmp(ftdata.label,'G2-MW',5);
yind = strncmp(ftdata.label,'G2-DJ',5);


X = ftdata.trial{1}(xind,:);
Y = ftdata.trial{1}(yind,:);
[X,Y] = remnan(X,Y);

figure
[c,lags] = crosscorr(X,Y,1000)
plot(1000*(lags./ftdata.fsample),c)
xlabel('Time lag (ms)'); ylabel('Cross-Correlation')