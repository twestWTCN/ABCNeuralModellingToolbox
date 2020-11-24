function [F meannpd] = constructGrangeMat(data,chloc_name,chlist,fsamp,N,R)
if isempty(N)
    N = floor(fsamp/R.obs.csd.df);
end
% N = R.obs.csd.pow2;
F_scale = R.frqz;
% Construct NPD matrix from data - take mean across channel replicates
for chloc = 1:size(chloc_name,2)
    chinds{chloc} = strmatch(chloc_name{chloc},chlist);
end
ftdata = [];
ftdata.label = R.chsim_name;
ftdata.trial = {data};
ftdata.fsample = fsamp;
ftdata.time{1} = linspace(0,numel(data(1,:))/fsamp,numel(data(1,:)));

cfg = [];
cfg.length = 0.25;
ftdata = ft_redefinetrial(cfg,ftdata);

cfg           = [];
cfg.foilim    = [min(F_scale) max(F_scale)];
cfg.method    = 'mtmfft';
cfg.taper     = 'hanning';
cfg.output    = 'fourier';
freq          = ft_freqanalysis(cfg, ftdata);

cfg           = [];
cfg.method    = 'instantaneous_causality';
instgranger           = ft_connectivityanalysis(cfg, freq);
cfg           = [];
cfg.method    = 'granger';
granger           = ft_connectivityanalysis(cfg, freq);
                    F = granger.freq;

for chI = 1:size(chloc_name,2)
    for chJ = 1:size(chloc_name,2)
        for p = 1:size(chinds{chI},1)
            chindsP = chinds{chI};
            for r = 1:size(chinds{chJ},1)
                chindsR = chinds{chJ};
                if chI == chJ
                    Pxy =   squeeze(mean(abs(freq.fourierspctrm(:,p,:)),1))';
                     Pxy = Pxy.*welchwin(length(Pxy))';
                    xcsd(p,r,1:3,:) = repmat(Pxy,3,1);
                else
                    
                    xcsd(p,r,1,:) = squeeze(instgranger.instantspctrm(chI,chJ,:)).*welchwin(length(F));
                    xcsd(p,r,2,:) = squeeze(granger.grangerspctrm(chI,chJ,:)).*welchwin(length(F));
                    xcsd(p,r,3,:) = squeeze(granger.grangerspctrm(chJ,chI,:)).*welchwin(length(F));
                    %                     zl = [10 11 12];
                    %                     for z = 1:3
                    %                         %                     [Pxy,F] = cpsd(data(chindsP(p),:),data(chindsR(r),:),hanning(N),[],N,fsamp);
                    %                         Pxy = f13(:,zl(z));
                    %                         if nargin>5
                    %                             Pxy = interp1(F,Pxy,F_scale);
                    %                         else
                    %                             Pxy =  Pxy(F>4);
                    %                         end
                    %                         Pxy = Pxy.*welchwin(length(Pxy))';
                    %                         xcsd(p,r,z,:) = Pxy;
                    %                     end
                end
            end
        end
        meannpd(chI,chJ,:,:) = (mean(mean(xcsd,1),2));
        clear xcsd
    end
end
% % take out symetrical CSD
% diaginds = [2 1; 3 1; 3 2];
% for i = 1:3
%     meancsd(diaginds(i,1),diaginds(i,2),:) = zeros(1,size(meancsd,3));
% end

end