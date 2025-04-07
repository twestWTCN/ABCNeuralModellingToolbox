function [r2mean,errorVec] = compareData_031121(R,sim_dat)

if ~isfield(R.objfx,'errorFx')
    R.objfx.errorFx = @fxSSE;
end

for empi = 1:numel(R.chdat_name)
    sim2Emp(empi) = find(strcmp(R.chsim_name,R.chdat_name{empi})); % these are the indices of the simdata
end
for dt = 1:numel(R.data.datatype)
    DatEmp  = R.data.feat_emp{dt}; % empirical
    DatSim  = sim_dat{dt}; % simulated
    
    switch R.data.datatype{dt}
        %% CSD
        case {'CSD','NPD'}
            for C = 1:numel(R.condnames)
                ic = 0; % these are the indices of the data
                for i = sim2Emp
                    ic = ic + 1;
                    jc = 0;
                    for j = sim2Emp
                        jc = jc + 1;
                        switch R.objfx.feattype
                            case 'complex'
                                if i~=j
                                    yfx = (squeeze(imag(DatSim(C,ic,jc,1,:)))); % model
                                    ffx = (squeeze(imag(DatEmp(C,ic,jc,1,:)))); % data
                                    r(1) = R.objfx.errorFx(yfx,ffx);
                                    
                                    yfx = (squeeze(real(DatSim(C,i,j,1,:)))); % model
                                    ffx = (squeeze(real(DatEmp(C,ic,jc,1,:)))); % data
                                    r(2) = R.objfx.errorFx(yfx,ffx);
                                    r2loop(C,ic,jc) = mean(r);
                                    
                                else
                                    yfx = squeeze(abs(DatSim(C,ic,jc,1,:))); % model
                                    ffx = squeeze(abs(DatEmp(C,ic,jc,1,:))); % data
                                    r(1) = R.objfx.errorFx(yfx,ffx);
                                    r2loop(C,ic,jc) = r(1);
                                end
                            case 'imaginary'
                                if i~=j
                                    yfx = (squeeze(imag(DatSim(C,ic,jc,1,:)))); % model
                                    ffx = (squeeze(imag(DatEmp(C,ic,jc,1,:)))); % data
                                    r(1) = R.objfx.errorFx(yfx,ffx);
                                    r2loop(C,ic,jc) = r(1); %mean(r);
                                    
                                else
                                    yfx = squeeze(abs(DatSim(C,ic,jc,1,:))); % model
                                    ffx = squeeze(abs(DatEmp(C,ic,jc,1,:))); % data
                                    r(1) = R.objfx.errorFx(yfx,ffx);
                                    r2loop(C,ii,jc) = r(1);
                                end
                            case 'absolute'
                                yfx = squeeze(abs(DatSim(C,ic,jc,1,:))); % model
                                ffx = squeeze(abs(DatEmp(C,ic,jc,1,:))); % data
                                r(1) = R.objfx.errorFx(yfx,ffx);
                                r2loop(C,ic,jc) = r(1);
                            case 'magnitude'
                                yfx = squeeze((DatSim(C,ic,jc,1,:))); % model
                                ffx = squeeze((DatEmp(C,ic,jc,1,:))); % data
                                r(1) = R.objfx.errorFx(yfx,ffx);
                                r2loop(C,ic,jc) = r(1);  %r(1);
                        end
                    end
                end
            end
            
            switch R.objfx.specspec
                case 'auto'
                    for C = 1:numel(R.condnames)
                        r2mean(C) = nanmean(diag(squeeze(r2loop(C,:,:))));
                    end
                    r2mean(dt) = mean(r2mean);
                    
                case 'cross'
                    for C = 1:numel(R.condnames)
                        r2C = squeeze(r2loop(C,:,:));
                        X = triu(r2C);
                        X(X==0) = [];
                        r2mean(C) = nanmean(X);
                    end
                    r2mean(dt) = mean(r2mean);
                case 'cross_only'
                    for C = 1:numel(R.condnames)
                        r2C = squeeze(r2loop(C,:,:));
                        r2mean(C) = nanmean(r2C(logical(~eye(j).*(triu(r2C)~=0))));
                    end
                    r2mean(dt) = nanmean(r2mean);
                case 'npd_only'
                    for C = 1:numel(R.condnames)
                        r2C = squeeze(r2loop(C,:,:));
                        r2mean(C) = nanmean(r2C(logical(~eye(j))));
                    end
                    r2mean(dt) = nanmean(r2mean);
                case 'npd'
                    for C = 1:numel(R.condnames)
                        r2C = squeeze(r2loop(C,:,:));
                        r2mean(C) = nanmean(r2C(:));
                    end
                    r2mean(dt) = nanmean(r2mean);
                    
            end
            
        case 'time' % time courses
                        TCemp  = R.data.feat_emp{dt}; % empirical
                        TCsim  = sim_dat{dt}; % simulated
            
                        for i = 1:size(TCemp,1)
                            try
                                yfx = TCsim(i,:);
                                ffx = TCemp(i,:);

                                if size(yfx,1)<size(yfx,2)
                                    yfx = yfx';
                                end
                                if size(ffx,1)<size(ffx,2)
                                    ffx = ffx';
                                end
                                r = R.objfx.errorFx(yfx,ffx);
                                r2loop(i) = r;
                            catch
                                r2loop(i) = -32;
                            end
                        end
                        r2mean(dt) = mean(r2loop);
        case 'none'
            r2mean(dt) = NaN;
        case {'FANO','DURPDF','BRSTPROF','ENVPDF','INTPDF'}
            r2loop = [];
            for C = 1:numel(R.condnames)
                ic = 0; % these are the indices of the data
                for i = sim2Emp
                    ic = ic + 1;
                    r2loop(:,C) = R.objfx.errorFx(DatSim(:,ic),DatEmp(:,ic));
                end
            end
            r2mean(dt) = nanmean(r2loop(:));
            %             fprintf('Fano error is: %0.3f  ',r2mean(dt))
    end
end

% Option to weight the respective features
r2mean(isnan(r2mean)) = -inf;
if isfield(R.objfx,'featweight')
    r2mean = r2mean.*R.objfx.featweight;
    r2mean(r2mean==0) = nan;
end
errorVec = [nanmean(r2mean(:)) r2mean]; %sum(r2mean);
r2mean = nanmean(r2mean(:));

