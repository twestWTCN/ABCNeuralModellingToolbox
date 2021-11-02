function a=  plotSweepSpectra(Hz,feat,featemp,cmap,legn,legsel,condsel,chsel,stcmap,scl)
if nargin<8
    chsel = 4;
end

% Mark out frequency borders
% p = patch([14 14 21 21],[1e-32 1e32 1e32 1e-32],'b');
% hold on
% p = patch([21 21 30 30],[1e-32 1e32 1e32 1e-32],'r');
% plot([14 1
% plot([14 14],[0 1],'k--')
% hold on
% plot([21 21],[0 1],'k--')
% plot([30 30],[0 1],'k--')

j = 0;
for i = condsel
    j = j+1;
    if max(squeeze(feat{i}{1}(1,4,4,1,:)));%<1e-4; % If STN is over reasonable level
        a(i) = plot(Hz,squeeze(feat{i}{1}(1,chsel(1),chsel(2),chsel(3),:))*scl,'color',cmap(i,:),'LineWidth',2);
        %     a(i) = plot(Hz,squeeze(feat{i}(1,4,1,2,:)),'color',cmap(i,:),'LineWidth',2);
        hold on
        pind(j) = i;
    end
end
a(i+1) = plot(Hz,squeeze(featemp{1}(1,chsel(1),chsel(2),chsel(3),:))*scl,'k','LineWidth',2);
a(i+2) = plot(Hz,squeeze(feat{2}{1}(1,chsel(1),chsel(2),chsel(3),:))*scl,'Color',stcmap(1,:),'LineWidth',2);
a(i+3) = plot(Hz,squeeze(feat{31}{1}(1,chsel(1),chsel(2),chsel(3),:))*scl,'Color',stcmap(2,:),'LineWidth',2);

% a(legsel(2)).LineStyle = '--';
% legend(a(legsel),legn)
xlim([4 38])
xlabel('Frequency (Hz)')
ylabel('Amplitude (uV Hz^-1)')
title('Simulated STN Spectra')
box off
grid on