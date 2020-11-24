function [F meannpd] = construct_spectral_ARGrangeMat(data,chloc_name,chlist,fsamp,N,R,normnoise)
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
ftdata.label = R.chsim_name ;
ftdata.trial = {data};
ftdata.fsample = fsamp;
ftdata.time{1} = linspace(0,numel(data(1,:))/fsamp,numel(data(1,:)));

cfg = [];
cfg.length = 0.5;
ftdata = ft_redefinetrial(cfg,ftdata);

cfg         = [];
cfg.order   = 12;
cfg.toolbox = 'bsmart';
mdata       = ft_mvaranalysis(cfg, ftdata);

cfg        = [];
cfg.method = 'mvar';
mfreq      = ft_freqanalysis(cfg, mdata);

cfg           = [];
cfg.method    = 'instantaneous_causality';
instgranger           = ft_connectivityanalysis(cfg, mfreq);
cfg           = [];
cfg.method    = 'granger';
granger           = ft_connectivityanalysis(cfg, mfreq);
                    F = granger.freq;
                    
% Compute Normalisation                    
    ftdata = [];
ftdata.label = R.chsim_name ;
ftdata.trial = {normnoise};
ftdata.fsample = fsamp;
ftdata.time{1} = linspace(0,numel(normnoise(1,:))/fsamp,numel(normnoise(1,:)));

cfg = [];
cfg.length = 0.5;
ftdata = ft_redefinetrial(cfg,ftdata);

cfg         = [];
cfg.order   = 12;
cfg.toolbox = 'bsmart';
mdata       = ft_mvaranalysis(cfg, ftdata);

cfg        = [];
cfg.method = 'mvar';
mfreqN      = ft_freqanalysis(cfg, mdata);

cfg           = [];
cfg.method    = 'instantaneous_causality';
instgrangerN           = ft_connectivityanalysis(cfg, mfreqN);
cfg           = [];
cfg.method    = 'granger';
grangerN           = ft_connectivityanalysis(cfg, mfreqN);
                    F = granger.freq;                
                    
chlen = size(chloc_name,2);

for chI = 1:size(chloc_name,2)
    for chJ = 1:size(chloc_name,2)
        for p = 1:size(chinds{chI},1)
            chindsP = chinds{chI};
            for r = 1:size(chinds{chJ},1)
                chindsR = chinds{chJ};
                if chI == chJ
                    Pxy =   squeeze(mean(abs(mfreq.crsspctrm(:,chI,:)),1))';
                     nPxy =   squeeze(mean(abs(mfreqN.crsspctrm(:,1,:)),1))';
%                      Pxy = Pxy.*welchwin(length(Pxy))';
                    Pxy = Pxy./max(nPxy);
                    xcsd(p,r,1:3,:) = repmat(Pxy,3,1);
                else
                    Cxy = squeeze(instgranger.instantspctrm(chI,chJ,:));
%                     nCxy =squeeze(instgrangerN.instantspctrm(chI,chJ,:));
%                     Cxy = (Cxy./max(nCxy));
                    xcsd(p,r,1,:) = Cxy; %.*welchwin(length(F));
                    
                    Cxy = squeeze(granger.grangerspctrm(chI,chJ,:));
                    nCxy = squeeze(grangerN.grangerspctrm(2,1,:));
                    Cxy = (Cxy./max(nCxy));
                    xcsd(p,r,2,:) = Cxy; %.*welchwin(length(F));
                    
                    Cxy = squeeze(granger.grangerspctrm(chJ,chI,:));
                    nCxy = squeeze(grangerN.grangerspctrm(2,1,:));
                    Cxy = (Cxy./max(nCxy));
                    xcsd(p,r,3,:) = Cxy;%.*welchwin(length(F));
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