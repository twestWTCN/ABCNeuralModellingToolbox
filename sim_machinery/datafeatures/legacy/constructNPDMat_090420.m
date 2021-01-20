function [F feat_out wflag meanconf] = constructNPDMat_090420(dataS,datinds,fsamp,N,R)
if isempty(N)
    N = floor(fsamp/R.obs.csd.df);
end

if ~isfield(R.obs,'logscale')
    R.obs.logscale = 1;
end

wflag = 0;
N = fix(N); %Ensure integer

% Construct NPD matrix from data - take mean across channel replicates
O = numel(R.condnames);
xcsd = nan(O,numel(R.siminds),numel(R.siminds),4,numel(R.frqz)); %initialize with empties
for C = 1:O
    for chI = datinds
        for chJ = datinds
            if chI == chJ
                [Pxy,F] = pwelch(squeeze(dataS{C}(chI,:)),hanning(2^N),[],2^(N),fsamp);
                F_scale = R.frqz;
                
                if nargin>4
                    Pxy = interp1(F,Pxy,F_scale,R.obs.trans.interptype);
                else
                    Pxy =  Pxy(F>4);
                end
                if R.obs.trans.logdetrend == 1
                    Pxy(Pxy<0) = 1e-32;
                    Pxy = log10(Pxy); F_scale = log10(F_scale);
                    [xCalc yCalc b Rsq] = linregress(F_scale',Pxy');
                    [dum bi] = intersect(F_scale,xCalc);
                    Pxy = Pxy(1,bi)-yCalc';
                    Pxy = 10.^Pxy; F_scale = 10.^(F_scale(bi));
                else
                    Pxy(isnan(F_scale)) = [];
                    F_scale(isnan(F_scale)) = [];
                end
                
                if R.obs.logscale == 1
                    Pxy = log10(Pxy);
                end
                if R.obs.trans.norm == 1
                    Pxy = (Pxy-nanmean(Pxy))./nanstd(Pxy);
                end
                
                if R.obs.trans.zerobase == 1
                    Pxy = Pxy - min(Pxy);
                end
                
                if R.obs.trans.gauss3 == 1
                    %                             Pxy = smoothdata(Pxy,'gaussian');
                    f = fit(F_scale',Pxy','gauss3');
                    Pxy = f(F_scale)';
                end
                
                if R.obs.trans.gausSm > 0
                    gwid = R.obs.trans.gausSm/diff(F_scale(1:2)); % 10 Hz smoothing
                    Pxy = smoothdata(Pxy,'gaussian',gwid);
                end
                Pxy(isnan(Pxy)) = 0;
                Pxy = Pxy; %.*tukeywin(length(Pxy),0.25)';
                xcsd(C,chJ,chI,1:4,:) = repmat(Pxy,4,1);
            else
                [f13,t,cl]=sp2a2_R2(squeeze(dataS{C}(chI,:))',squeeze(dataS{C}(chJ,:))',fsamp,N);
                F = f13(:,1);
                zl = [10 11 12 13];
                for z = 1:4
                    %                     [Pxy,F] = cpsd(data(chindsP(p),:),data(chindsR(r),:),hanning(N),[],N,fsamp);
                    Pxy = f13(:,zl(z));
                    
                    F_scale = R.frqz;
                    F_scale(isnan(F_scale)) = [];
                    if nargin>4
                        Pxy = interp1(F,Pxy,F_scale,'pchip');
                    else
                        Pxy =  Pxy(F>4);
                    end
                    
                    if R.obs.trans.norm == 1
                        Pxy = (Pxy-nanmean(Pxy))./nanstd(Pxy);
                        Pxy = Pxy - min(Pxy);
                    end
                    
                    if R.obs.trans.gauss3 == 1
                        %                             Pxy = smoothdata(Pxy,'gaussian');
                        f = fit(F_scale',Pxy','gauss3');
                        Pxy = f(F_scale)';
                    end
                    if R.obs.trans.gausSm > 0
                        gwid = R.obs.trans.gausSm/diff(F_scale(1:2)); % 10 Hz smoothing
                        Pxy = smoothdata(Pxy,'lowess',gwid);
                    end
                    xcsd(chI,chJ,z,:) = Pxy;
                end
            end
        end
    end
end

feat_out = xcsd;
F = F_scale;
meanconf = [];