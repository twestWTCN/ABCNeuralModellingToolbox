function [F meannpd] = construct_spectral_NPGrangeMat(data,chloc_name,chlist,fsamp,N,R,normnoise)
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
ftdata.label = R.chsim_name; % 'N1' 'N2'];
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


% Normalization
ftdata = [];
ftdata.label = R.chsim_name; % 'N1' 'N2'];
ftdata.trial = {normnoise};
ftdata.fsample = fsamp;
ftdata.time{1} = linspace(0,numel(normnoise(1,:))/fsamp,numel(normnoise(1,:)));

cfg = [];
cfg.length = 0.25;
ftdata = ft_redefinetrial(cfg,ftdata);

cfg           = [];
cfg.foilim    = [min(F_scale) max(F_scale)];
cfg.method    = 'mtmfft';
cfg.taper     = 'hanning';
cfg.output    = 'fourier';
freqN          = ft_freqanalysis(cfg, ftdata);

cfg           = [];
cfg.method    = 'instantaneous_causality';
instgrangerN   = ft_connectivityanalysis(cfg, freqN);
cfg           = [];
cfg.method    = 'granger';
grangerN       = ft_connectivityanalysis(cfg, freqN);

F = granger.freq;
chlen = size(chloc_name,2);
for chI = 1:size(chloc_name,2)
    for chJ = 1:size(chloc_name,2)
        for p = 1:size(chinds{chI},1)
            chindsP = chinds{chI};
            for r = 1:size(chinds{chJ},1)
                chindsR = chinds{chJ};
                if chI == chJ
                    Pxy =   squeeze(mean(abs(freq.fourierspctrm(:,p,:)),1))';
                    nPxy =   squeeze(mean(abs(freqN.fourierspctrm(:,1,:)),1))';
%                     Pxy = Pxy.*welchwin(length(Pxy))';
                    Pxy = Pxy./max(nPxy);
                    xcsd(p,r,1:3,:) = repmat(Pxy,3,1);
                else
                    Cxy = squeeze(instgranger.instantspctrm(chI,chJ,:));
                    %                     nCxy =squeeze(instgranger.instantspctrm(chI+chlen,chJ+chlen,:));
                    %                     Cxy = (Cxy./mean(nCxy));
                    xcsd(p,r,1,:) = Cxy;%.*welchwin(length(F));
                    
                    Cxy = squeeze(granger.grangerspctrm(chI,chJ,:));
                    nCxy = squeeze(grangerN.grangerspctrm(2,1,:));
                    Cxy = (Cxy./max(nCxy));
                    xcsd(p,r,2,:) = Cxy;%.*welchwin(length(F));
                    
                    Cxy = squeeze(granger.grangerspctrm(chJ,chI,:));
                    nCxy = squeeze(grangerN.grangerspctrm(2,1,:)); % Always same direction as no reverse
                    Cxy = (Cxy./max(nCxy));
                    xcsd(p,r,3,:) = Cxy; %.*welchwin(length(F));
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