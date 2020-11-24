function plotABCSpectraOnly(Hz,EMP,SIM)
for i = 1:size(SIM,2)
    if i<5 && ~isempty(EMP)
        emp = squeeze(EMP(1,i,i,1,:));
    else
        emp = nan;
    end
    
    sim = squeeze(SIM(1,i,i,1,:));
    
    subplot(1,6,i)
    x = Hz;
    plot(x,sim)
    hold on
    plot(x,emp)
    axis square
    xlabel('Frequency (Hz)')
    ylabel('Normalized Power')
    ylim([0 6])
end