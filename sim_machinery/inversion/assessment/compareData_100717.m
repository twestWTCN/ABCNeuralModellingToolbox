function r2mean = compareData_100717(R,sim_dat)
switch R.data.datatype
    %% CSD
    case 'CSD'
        NPDemp  = R.data.feat_emp; % empirical
        NPDsim  = sim_dat; % simulated
        for C = 1:numel(R.condnames)
            for i = 1:size(NPDemp,2)
                for j = 1:size(NPDemp,3)
                    switch R.objfx.feattype
                        case 'complex'
                            if i~=j
                                yfx = (squeeze(imag(NPDsim(C,i,j,1,:))));
                                ffx = (squeeze(imag(NPDemp(C,i,j,1,:))));
                                r(1) = -RMSE(yfx,ffx);
                                
                                yfx = (squeeze(real(NPDsim(C,i,j,1,:))));
                                ffx = (squeeze(real(NPDemp(C,i,j,1,:))));
                                r(2) = -RMSE(yfx,ffx);
                                r2loop(C,i,j) = mean(r);
                                
                            else
                                yfx = squeeze(abs(NPDsim(C,i,j,1,:)));
                                ffx = squeeze(abs(NPDemp(C,i,j,1,:)));
                                r(2) = -RMSE(yfx,ffx);
                                
                                r2loop(C,i,j) = r(2);
                            end
                        case 'imaginary'
                            if i~=j
                                yfx = (squeeze(imag(NPDsim(C,i,j,1,:))));
                                ffx = (squeeze(imag(NPDemp(C,i,j,1,:))));
                                r(1) = goodnessOfFit(yfx,ffx,'NRMSE');
                                
                                r2loop(C,i,j) = r(1); %mean(r);
                                
                            else
                                yfx = squeeze(abs(NPDsim(C,i,j,1,:)));
                                ffx = squeeze(abs(NPDemp(C,i,j,1,:)));
                                r(2) = goodnessOfFit(yfx,ffx,'NRMSE');
                                
                                r2loop(C,i,j) = r(2);
                            end
                        case 'absolute'
                            yfx = squeeze(abs(NPDsim(C,i,j,1,:)));
                            ffx = squeeze(abs(NPDemp(C,i,j,1,:)));
                            r(1) = goodnessOfFit(yfx,ffx,'NRMSE');
                            r2loop(C,i,j) = r(1);
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
                r2mean = median(r2mean);
                
            case 'cross'
                for C = 1:numel(R.condnames)
                    r2C = squeeze(r2loop(C,:,:));
                    X = triu(r2C);
                    X(X==0) = [];
                    %                     r2mean(C) = nanmean(r2C(triu(r2C)~=0));
                    r2mean(C) = sum(X);
                end
                %         r2mean = sum(r2loop(triu(r2loop)~=0));
                r2mean = sum(r2mean);
                %                 simdat = yfxx(:); simdat(isnan(simdat)) = 0;
                %                 empdat = ffxx(:); empdat(isnan(empdat)) = 0;
                %                 r2mean = goodnessOfFit(simdat,empdat,'NRMSE');
            case 'cross_only'
                for C = 1:numel(R.condnames)
                    r2C = squeeze(r2loop(C,:,:));
                    r2mean(C) = nanmedian(r2C(logical(~eye(j).*(triu(r2C)~=0))));
                end
                r2mean = median(r2mean);
        end
        %% NPD
    case 'NPD'
        NPDemp  = R.data.feat_emp; % empirical
        NPDsim  = sim_dat; % simulated
        yfxx = [];
        ffxx = [];
        r = [];
        for C = 1:numel(R.condnames)
            for i = 1:size(NPDemp,2)
                for j = 1:size(NPDemp,2)
                    switch R.objfx.feattype
                        case 'ForRev'
                            if i==j
                                yfx = (squeeze(NPDsim(C,i,j,1,:)));
                                ffx = (squeeze(NPDemp(C,i,j,1,:)));
                                ffx(isnan(yfx)) = [];yfx(isnan(yfx)) = [];
                                yfx(isnan(ffx)) = [];ffx(isnan(ffx)) = [];
                                
                                if size(yfx,1)>size(yfx,2);yfx = yfx'; end
                                if size(ffx,1)>size(ffx,2);ffx = ffx'; end
                                yfxx = [yfxx yfx];
                                ffxx  = [ffxx ffx];
                                yfx = yfx+1; ffx = ffx+1;
                                %                                 r(2) = goodnessOfFit(yfx',ffx','NRMSE');
                                %                                 r(2) =  1-(RMSE(yfx',ffx')./std(ffx));
                                r(2) = rsquare(yfx',ffx');
                                r2loop(C,i,j) = r(2);
                            elseif j>i
                                for k = 2:3 % 1 is zerolag
                                    yfx = (squeeze(NPDsim(C,i,j,k,:)));
                                    ffx = (squeeze(NPDemp(C,i,j,k,:)));
                                    
                                    if size(yfx,1)>size(yfx,2);yfx = yfx'; end
                                    if size(ffx,1)>size(ffx,2);ffx = ffx'; end
                                    
                                    ffx(isnan(yfx)) = [];yfx(isnan(yfx)) = [];
                                    yfx(isnan(ffx)) = [];ffx(isnan(ffx)) = [];
                                    
                                    yfxx = [yfxx yfx];
                                    %                                     ffxx  = [ffxx ffx];
                                    yfx = yfx+1; ffx = ffx+1;
                                    %                                     r(k) = goodnessOfFit(yfx',ffx','NRMSE');
                                    %                                     r(k) =  1-(RMSE(yfx',ffx')./std(ffx));
                                    r(k) = rsquare(yfx',ffx');
                                end
                                r2loop(C,i,j) = mean(r(2));
                                r2loop(C,j,i) = mean(r(3));
                            end
                    end
                end
            end
        end
        switch R.objfx.specspec
            case 'auto'
                for C = 1:numel(R.condnames)
                    r2mean(C) = nanmedian(diag(squeeze(r2loop(C,:,:))));
                end
                r2mean = median(r2mean);
                
            case 'cross'
                for C = 1:numel(R.condnames)
                    r2C = squeeze(r2loop(C,:,:));
                    %                     r2mean(C) = nanmean(r2C(triu(r2C)~=0));
                    r2mean(C) = nanmedian(r2C(:));
                end
                %         r2mean = sum(r2loop(triu(r2loop)~=0));
                r2mean = median(r2mean);
                %                 simdat = yfxx(:); simdat(isnan(simdat)) = 0;
                %                 empdat = ffxx(:); empdat(isnan(empdat)) = 0;
                %                 r2mean = goodnessOfFit(simdat,empdat,'NRMSE');
            case 'cross_only'
                for C = 1:numel(R.condnames)
                    r2C = squeeze(r2loop(C,:,:));
                    r2mean(C) = nanmedian(r2C(logical(~eye(j).*(triu(r2C)~=0))));
                end
                r2mean = median(r2mean);
        end
        %% TIME
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
        r2mean = mean(r2loop);
    case 'none'
        r2mean = NaN;
        
end