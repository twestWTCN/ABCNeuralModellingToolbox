function plotTimeSeriesGen(data,fsample,legname,condname)
i = 0;
step = 6;
clf
for cond = 1:numel(condname)
    for rep = 1:2
        i = i+1;
        subplot(2,2,i)
        for ch = 1:size(data{cond},1)
            X = data{cond}(ch,:);
            X = normalize(X);
            X = X-(step*(ch-1));
            t = linspace(0,size(X,2)./fsample,size(X,2));
            plot(t,X)
            hold on
        end
        if rep ==2
            xlim([t(fix(0.5*numel(t))) t(fix(0.55*numel(t)))])
        end
        xlabel('Time (s)')
        ylabel('Normalized Amplitude')
        title(condname{cond})
            grid on
            ylim([-25 10])
    end
    if i == 4
        legend(legname)
    end
end
