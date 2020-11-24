function a = plotObjectiveViolin(R,compStruc,shortlab)
violin(compStruc.r2repSave,'facecolor',R.plot.cmap,'medc','k:','xlabel',shortlab);
hold on
% plot([0 numel(R.modcomp.modN)+1],[R.modcomp.modEvi.epspop R.modcomp.modEvi.epspop],'k--')
xlabel('Model'); ylabel('NMRSE'); grid on; ylim([-2 1])
a = gca; a.XTick = 1:numel(R.modcomp.modN);
a.XTickLabel = shortlab;
