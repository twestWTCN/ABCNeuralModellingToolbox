function [xsims R] = observe_data_170118(xstore,m,p,R)

xsims = xstore(R.obs.states,:);
% Delete burnin
xsims(:,1:round(R.obs.brn*(1/R.IntP.dt))) = [];
tvec_obs = R.IntP.tvec;
tvec_obs(:,1:round(R.obs.brn*(1/R.IntP.dt))) = [];
R.IntP.tvec_obs = tvec_obs;

for i = 1:length(R.obs.gainmeth)
    switch R.obs.gainmeth{i}
        case 'leadfield'
            LF = R.obs.LF.*exp(p.obs.LF);
            LFF = zeros(m.m);
            LFF(eye(size(LFF))~=0) = LF;
            xsims = LFF*xsims;
        case 'unitvar'
            for i = 1:size(xsims,1)
                xsims(i,:) = (xsims(i,:) - mean(xsims(i,:)))./std(xsims(i,:));
            end
        case 'mixing'
            %% REPLACE WITH DISTANCE MATRIX
                mixdeg = R.obs.mixing(1).*exp(p.obs.mixing(1));
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
    end
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
