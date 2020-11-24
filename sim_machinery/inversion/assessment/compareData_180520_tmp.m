function r2mean = compareData_180520(R,sim_dat)
for empi = 1:numel(R.chdat_name)
    sim2Emp(empi) = find(strcmp(R.chsim_name,R.chdat_name{empi}));
end
for dt = 1:numel(R.data.datatype)
                NPDemp  = R.data.feat_emp{dt}; % empirical
            NPDsim  = sim_dat{dt}; % simulated

    switch R.data.datatype{dt}        
        %% CSD
        case {'CSD','NPD'}
            for C = 1:numel(R.condnames)
                ic = 0;
                for i = sim2Emp
                    ic = ic + 1;
                    jc = 0;
                    for j = sim2Emp
                        jc = jc + 1;
                        switch R.objfx.feattype
                            case 'complex'
                                if i~=j
                                    yfx = (squeeze(imag(NPDsim(C,i,j,1,:))));
                                    ffx = (squeeze(imag(NPDemp(C,i,j,1,:))));
                                    r(1) = -RMSE(yfx,ffx);
                                    
                                    yfx = (squeeze(real(NPDsim(C,i,j,1,:))));
                                    ffx = (squeeze(real(NPDemp(C,i,j,1,:))));
                                    r(2) = -RMSE(yfx,ffx);
                                    r2loop(C,ic,jc) = mean(r);
                                    
                                else
                                    yfx = squeeze(abs(NPDsim(C,i,j,1,:)));
                                    ffx = squeeze(abs(NPDemp(C,i,j,1,:)));
                                    r(2) = -RMSE(yfx,ffx);
                                    
                                    r2loop(C,ic,jc) = r(2);
                                end
                            case 'imaginary'
                                if i~=j
                                    yfx = (squeeze(imag(NPDsim(C,i,j,1,:))));
                                    ffx = (squeeze(imag(NPDemp(C,i,j,1,:))));
                                    r(1) = -RMSE(yfx,ffx);
                                    
                                    r2loop(C,ic,jc) = r(1); %mean(r);
                                    
                                else
                                    yfx = squeeze(abs(NPDsim(C,i,j,1,:)));
                                    ffx = squeeze(abs(NPDemp(C,i,j,1,:)));
                                    r(1) = -RMSE(yfx,ffx);
                                    
                                    r2loop(C,ii,jc) = r(1);
                                end
                            case 'absolute'
                                yfx = squeeze(abs(NPDsim(C,i,j,1,:)));
                                ffx = squeeze(abs(NPDemp(C,i,j,1,:)));
                                r(1) = -RMSE(yfx,ffx);
                                r2loop(C,ic,jc) = r(1);
                            case 'magnitude'
                                yfx = squeeze((NPDsim(C,i,j,1,:)));
                                ffx = squeeze((NPDemp(C,i,j,1,:)));
                                r(1) = -RMSE(yfx,ffx);
                                r2loop(C,ic,jc) = r(1);
                        end
                    end
                end
            end
            % r2loop = triu(r2loop);
            % r2loop = diag(r2loop);
            % r2loop = r2loop(1,1);
            % r2loop = reshape(r2loop,1,[]);
            % r2loop(isnan(r2loop)) = [];
            % r2mean = mean(r2loop,2);
            switch R.objfx.specspec
                case 'auto'
                    for C = 1:numel(R.nanmedian)
                        r2mean(C) = nanmean(diag(squeeze(r2loop(C,:,:))));
                    end
                    r2mean(dt) = mean(r2mean);
                    
                case 'cross'
                    for C = 1:numel(R.condnames)
                        r2C = squeeze(r2loop(C,:,:));
                        X = triu(r2C);
                        X(X==0) = [];
                        %                     r2mean(C) = nanmean(r2C(triu(r2C)~=0));
                        r2mean(C) = nanmean(X);
                    end
                    %         r2mean = sum(r2loop(triu(r2loop)~=0));
                    r2mean(dt) = mean(r2mean);
                    %                 simdat = yfxx(:); simdat(isnan(simdat)) = 0;
                    %                 empdat = ffxx(:); empdat(isnan(empdat)) = 0;
                    %                 r2mean = goodnessOfFit(simdat,empdat,'NRMSE');
                case 'cross_only'
                    for C = 1:numel(R.condnames)
                        r2C = squeeze(r2loop(C,:,:));
                        r2mean(C) = nansum(r2C(logical(~eye(j).*(triu(r2C)~=0))));
                    end
                    r2mean(dt) = nanmean(r2mean);
            end
            
        case 'time' % time courses
            TCemp  = R.data.feat_emp{1}; % empirical
            TCsim  = sim_dat{1}; % simulated
            
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
                    r = rsquare(yfx,ffx);
                    %  r = goodnessOfFit(yfx,ffx,'NRMSE');
                    r2loop(i) = r;
                catch
                    r2loop(i) = -32;
                end
            end
            r2mean(dt) = mean(r2loop);
        case 'none'
            r2mean(dt) = NaN;
        case {'FANO','DUR'}
            r2loop = [];
            for C = 1:numel(R.condnames)
                r2loop(:,C) = -RMSE(NPDemp,NPDsim(:,R.datinds));
            end
            r2mean(dt) = nanmean(r2loop)*50000;
            fprintf('Fano error is: %0.3f  ',r2mean(dt))
    end
end
r2mean = mean(r2mean(:)); %sum(r2mean);

