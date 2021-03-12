function [F feat_out wflag meanconf] = constructGenCrossMatrix(dataS,datinds,fsamp,N,R)
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
    
    switch R.data.datatype{1}
        case 'CSD'
            dataX = dataS{C}(datinds,:)';
            
            % Compute the CrossSpectral Density for Everything
            [csdMaster,fMaster] = cpsd(dataX,dataX,hanning(2^N),[],2^N,fsamp,'mimo');
            if numel(size(csdMaster))<3
                X = nan(size(csdMaster,1),2,2);
                X(:,1,1) = csdMaster;
                csdMaster = X;
            end
        case 'NPD'
            i = 0;
            for chI = 1:numel(datinds)
                i = i + 1;
                j = 0;
                for chJ = 1:numel(datinds)
                    j = j + 1;
                    [f13,t,cl]=sp2a2_R2(squeeze(dataS{C}(datinds(chJ),:))',squeeze(dataS{C}(datinds(chI),:))',fsamp,N);
                    fMaster = f13(:,1);
                    if chI == chJ
                        csdMaster(:,j,i) = 10.^(f13(:,2));
                    else
                        csdMaster(:,j,i) = f13(:,11)*R.obs.trans.npdscalar;
                    end
                end
            end
            
    end
    
    for i = 1:numel(datinds)
        for j = 1:numel(datinds)
            if i == j
                % Your Univariate Measure
                Pxy = squeeze(csdMaster(:,j,i));
                F_scale = fMaster;
                
                if R.obs.trans.logdetrend == 1
                    Pxy(Pxy<0) = 1e-32;
                    tailinds = (F_scale>68);
                    Pxy = log10(Pxy); F_scale = log10(F_scale);
                    [dum1 dum2 b Rsq] = linregress(F_scale(tailinds),Pxy(tailinds));
                    yCalc = [ones(length(F_scale),1) F_scale]*b;
                    Pxy = Pxy-yCalc;
                    Pxy = 10.^Pxy; F_scale = 10.^(F_scale);
                else
                    Pxy(isnan(F_scale)) = [];
                    F_scale(isnan(F_scale)) = [];
                end
                
                if R.obs.logscale == 1
                    Pxy = log10(Pxy);
                end
                
                
                if R.obs.trans.gauss3 == 1
                    %                             Pxy = smoothdata(Pxy,'gaussian');
                    f = fit(F_scale',Pxy','gauss3');
                    Pxy = f(F_scale)';
                end
                
                if nargin>4
                    Pxy = interp1(F_scale,Pxy,R.frqz,R.obs.trans.interptype);
                    Pxy = Pxy;
                else
                    Pxy =  Pxy(F>4);
                end
                
                if R.obs.trans.gausSm > 0
                    gwid = R.obs.trans.gausSm/diff(R.frqz(1:2)); % 10 Hz smoothing
                    Pxy = smoothdata(Pxy,'gaussian',gwid);
                end
                
                if R.obs.trans.norm == 1
                    Pxy = (Pxy-nanmean(Pxy))./nanstd(Pxy);
                end
                if R.obs.trans.zerobase == 1
                    Pxy = Pxy - min(Pxy);
                end
                
            elseif i ~= j % Diagonal
                % Your Functional Connectivity Metric
                Pxy = squeeze(csdMaster(:,j,i));
                F_scale = fMaster;
                
                if R.obs.trans.norm == 1
                    Pxy = (Pxy-nanmean(Pxy))./nanstd(Pxy);
                    Pxy = Pxy - min(Pxy);
                end
                
                if R.obs.trans.gauss3 == 1
                    %                             Pxy = smoothdata(Pxy,'gaussian');
                    f = fit(F_scale',Pxy','gauss3');
                    Pxy = f(F_scale)';
                end
                
                if nargin>4
                    Pxy = interp1(F_scale,Pxy,R.frqz,R.obs.trans.interptype);
                else
                    Pxy =  Pxy(F>4);
                end
                
                if R.obs.trans.gausSm > 0
                    gwid = R.obs.trans.gausSm/diff(R.frqz(1:2)); % 10 Hz smoothing
                    Pxy = smoothdata(Pxy,'gaussian',gwid);
                end
            end
            xcsd(C,j,i,1:4,:) = repmat(Pxy,4,1);
        end
    end
end
if R.obs.trans.normcat == 1
    % Normalize each component by concatanating the conditions
    for chI = 1:numel(datinds)
        for chJ = 1:numel(datinds)
            Xd = xcsd(:,chJ,chI,1:4,:); %here you select both conditions
            XM = mean(Xd(:));
            XV = std(Xd(:));
            Xd = xcsd(:,chJ,chI,1:4,:); %here you select both conditions
            xcsd(:,chJ,chI,1:4,:) = (Xd)./XV; % Rescale
        end
    end
end
feat_out{1} = xcsd;
F{1} = R.frqz;

for C = 1:O
    if numel(R.data.datatype)>1
        for fcnt = 2:numel(R.data.datatype)
            nbs = [];
            switch R.data.datatype{fcnt}
                case 'FANO'
                    dataX = dataS{C}(datinds,:)';
                    nb = []; E = [];
                    for i = 1:size(dataX,2)
                        X = bandpass(normalize(dataX(:,i))',[24 52],fsamp);
                        XH = abs(hilbert(X));
                        [nb(:,i),E(:,i)] = histcounts(XH,R.data.feat_xscale{fcnt},'Normalization','pdf');
                    end
                    
                case 'DURPDF'
                    dataX = dataS{C}(datinds,:)';
                    nb = []; E = [];
                    for i = 1:size(dataX,2)
                        [nb(:,i,C),nbs(:,i,C),E(:,i,C)] = burstDurHist(dataX(:,i)',[15 25],fsamp,R.data.feat_xscale{fcnt});
                    end
                case {'BRSTPROF','ENVPDF'}
                    dataX = dataS{C}(datinds,:)';
                    dataX = bandpass(dataX,[15 25],fsamp);
                    dataX = (dataX-mean(dataX))./std(dataX);
                    nAvg                = 1; % the time series is divided into nAvg segments to plot sem error bars
                    minBurstDuration    = 0.1; % burst are only considered if longer than this duration (in s)
                    xPerc               = R.data.feat_xscale{fcnt}; % vector of thresholds
                    
                    E = xPerc;
                    try
                        E = R.data.feat_xscale{strcmp(R.data.datatype,'BRSTPROF')};
                    catch
                        E = 5:85; % default percentiles
                    end
                    try
                        XPDF = R.data.feat_xscale{strcmp(R.data.datatype,'ENVPDF')};
                    catch
                        XPDF = linspace(0,max(env(:),100)); % default percentiles
                    end
                    
                    
                    for i = 1:size(dataX,2)
                        env = abs(hilbert(dataX(:,1)));
                        [avgBurstDuration,semBurstDuration,avgPDF,semPDF] = burstDurWrapper(env,E,nAvg,1/fsamp,minBurstDuration,[],XPDF); % burst profile
                    end
                    avgBurstDuration(isnan(nbs)) = -20; % IS THIS OK?
                    switch R.data.datatype{fcnt}
                        case 'BRSTPROF'
                            nbs(:,i,C) = avgBurstDuration;
                        case 'ENVPDF'
                            nbs(:,i,C) = avgPDF;
                    end
                    %                 nbs(isnan(nbs)) = -20; % IS THIS OK?
            end
            feat_out{fcnt} = nbs;
            F{fcnt} = E;
        end
    end
end

meanconf = [];
