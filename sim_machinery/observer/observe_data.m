function [xsims_c R wflag] = observe_data(xstore,m,p,R)
wflag = 0;

for condsel = 1:numel(R.condnames)
    xsims = xstore{condsel}(R.obs.outstates,:); % select simulation states
    %     xsims = xsims(R.obs.obsstates,:); % and states to be observed (or put into LF)
    
    % Delete burnin
    if size(xsims,2) > 5*round(R.obs.brn*(1/R.IntP.dt))
        xsims(:,1:round(R.obs.brn*(1/R.IntP.dt))) = [];
    else
        wflag = 1;
        xsims_c{condsel} = xsims;
        warning('Simulation is shorter than burn length!!!')
        return
    end
    if any(isnan(xsims(:)))
        xsims_c{condsel} = xsims;
        warning('Simulation contains NaNs!!!')
        return
    end
    tvec_obs = R.IntP.tvec;
    tvec_obs(:,1:round(R.obs.brn*(1/R.IntP.dt))) = [];
    R.IntP.tvec_obs = tvec_obs;
    
    for i = 1:length(R.obs.gainmeth)
        switch R.obs.gainmeth{i}
            case 'obsnoise'
                CN = (R.obs.Cnoise.*exp(p.obs.Cnoise))';
                xsims = xsims + randn(size(xsims)).*CN.*std(xsims,[],2); % SNR of target
            case 'obsCnoise'
                CN = (R.obs.Cnoise.*exp(p.obs.Cnoise))';
                alpha = R.obs.AlpNoise.*exp(p.obs.AlpNoise)';
                for ij = 1:size(xsims,1)
                    U = ffGn(size(xsims,2),(alpha+1)/2, std(xsims(ij,:)), 0).*CN;
                    xsims(ij,:) = xsims(ij,:) + U;
                end
            case 'leadfield'
                LF = m.obs.LF.*exp(p.obs.LF);
                sigmix = repmat(1-LF,m.m,1).*eye(m.n);
                sigmix = sigmix + (repmat(LF/(m.n-1),m.n,1).*~eye(m.n));
                xsims = sigmix*xsims;
            case 'difference'
                xsims = [xsims(:,1) diff(xsims,1,2)];
            case 'unitvar'
                for j = 1:size(xsims,1)
                    xsims(j,:) = (xsims(j,:) - mean(xsims(j,:)))./std(xsims(j,:));
                end
            case 'unitvarConcat'
                LM = [];
                for C = 1:numel(xstore)
                    LM = [LM xstore{C}(R.obs.outstates,round(R.obs.brn*(1/R.IntP.dt)):end)];
                end
                XM = mean(LM,2);
                XV = std(LM,[],2);
                xsims = (xsims-XM)./XV;
            case 'mixing'
                %% REPLACE WITH DISTANCE MATRIX
                mixdeg = m.obs.mixing(1).*exp(p.obs.mixing(1));
                sigmix = repmat(1-mixdeg,m.m,1).*eye(m.m);
                sigmix = sigmix + (repmat(mixdeg/(m.m-1),m.m,1).*~eye(m.m));
                xsims = sigmix*xsims;
            case 'submixing'
                %% REPLACE WITH DISTANCE MATRIX
                mixdeg = R.obs.mixing.*exp(p.obs.mixing);
                cortmix = [mixdeg(1) repmat(mixdeg(2),1,m.m-1)];
                cortmix = [1-(mixdeg(1)*(m.m-1)) repmat(mixdeg(1),1,m.m-1)];
                submix = ~eye(m.m-1,m.m).*repmat(mixdeg(1),m.m-1,m.m);
                submix(logical(eye(size(submix)))) = 1-(mixdeg(2));
                
                mix = [cortmix; circshift(submix,1,2)];
                %             m.m = 6;
                %             dm = eye(m.m) + (~eye(m.m).*repmat(mixdeg(2),m.m,m.m));
                %             dm(2:m.m,1) = repmat(mixdeg(1),1,m.m-1);
                %             dm(1,2:m.m) = repmat(mixdeg(1),m.m-1,1);
                %             dm(logical(eye(size(dm)))) = 0;
                %             dm = dm.*0.001
                xsims = mix*xsims;
            case 'lowpass'
                for i = 1:size(xsims,1)
                    x = xsims(i,:);
                    xsims(i,:) = filtfilt(R.obs.lowpass.fwts,1,x);
                end
            case 'highpass'
                for i = 1:size(xsims,1)
                    x = xsims(i,:);
                    xsims(i,:) = filtfilt(R.obs.highpass.fwts,1,x);
                end
            case 'boring'
                %                 figure(100);
                %                 clf
                %                 plot(xsims'); shg
                %
                montoncheck = [];
                for j = 1:size(xsims,1)
                    swX = slideWindow(xsims(j,:), floor(size(xsims(j,:),2)/3), 0);
                    for swin = 1:size(swX,2)
                        A = swX(:,swin);
                        Env = abs(hilbert(A));
                        [acf,lags] = xcorr(Env,Env,1000,'coeff');
                        %                         fft(acf)
                        acfEnvcheck(j,swin) = any(abs(acf(abs(lags)>100))>0.85);
                        [acf,lags] = xcorr(acf,acf,1000,'coeff');
                        acf2Envcheck(j,swin) = any(abs(acf(abs(lags)>300))>0.95);
                        
                    end
                    Env = abs(hilbert(xsims(j,:)));
                    % Subsample
                    Env = Env(1:100:end);
                    [tau pc] = corr((1:size(Env,2))',Env','type','Kendall');
                    montoncheck(j) = pc<0.05;
                    %                     Xstab(i) = std(diff(abs(hilbert(xsims(i,:)))))<0.005;
                end
                
                if any(acfEnvcheck(:)) || any(acf2Envcheck(:))
                    disp('SimFx is perfectly periodic!!')
                    wflag = 1;
                elseif sum(montoncheck(:))>2 %any(acfcheck) ||
                    disp('SimFx seems unstable!')
                    wflag = 1;
                end
                if wflag == 0
                    a = 1;
                end
                %                 pause(2)
            case 'FANO'
                
                R.sim.fano(:,condsel) = computeFano(xsims,1/R.IntP.dt);
            case 'flatTail'
                
                [f,fx] =  pwelch(xsims',1/R.IntP.dt,[],1/R.IntP.dt,1/R.IntP.dt);
                tailSNR =  10.*log10(sum(f(fx<=48,:))./sum(f(fx>48 & fx<248,:)));
                if any(tailSNR< 1)
                    wflag = 1;
                    disp('TailSNR failed')
                end
            case 'symmetrical'
                for j = 1:size(xsims,1)
                    X = (xsims(j,:) - mean(xsims(j,:)))./std(xsims(j,:));
                    XL = abs(X(X<0)); XR = abs(X(X>0));
                    skew(j) = abs(mean(XL)-mean(XR));
                end
%                 
%                 figure(100);
%                 clf
%                 plot(xsims');
%                 xlim([1 1.1]*10^4)
%                 shg
                
                if any(skew>0.5)
                    wflag = 1;
                    disp('Siganals not symmetrical')
                end
                
        end
    end
    xsims_c{condsel} = xsims;
end

% if R.obs.condchecker == 1
% %     xsims_c{1} = xsims_c{1} - mean(xsims_c{1},2);
% %         xsims_c{2} = xsims_c{2} - mean(xsims_c{2},2);
%     maxD = ((max(xsims_c{1},[],2)-max(xsims_c{2},[],2))./max(xsims_c{2},[],2)).*100;
%        if any(abs(maxD)>1e4)
%            wflag = 1;
%                   disp('Scale Differeence of conditions is too large!')
%
%        end
% end

if wflag == 0
    a = 1;
    %                     plot([xsims_c{1} xsims_c{2}]')
end
% if R.obs.norm
%     % Normalise
%     for i = 1:size(xstore,1)
%         xtmp = xstore(i,:);
%         xtmp = (xtmp-mean(xtmp))/std(xtmp);
%         %         xtmp = xtmp.*(m(i).lfpgain*exp(p(i).lfpgain)); % Multiply by gainfield
%         xsims(i,:) = xtmp;
%     end
% end

% Amplify to gain
