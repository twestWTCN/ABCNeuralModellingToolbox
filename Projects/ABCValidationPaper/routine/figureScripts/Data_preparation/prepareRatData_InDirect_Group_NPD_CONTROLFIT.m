function R = prepareRatData_InDirect_Group_NPD_CONTROLFIT(R,plotop,nulldat)
if nargin<2
    plotop = 1;
end
if nargin<3
    nulldat = 0;
end
load([R.path.datapath '\NPD_paper_RatNPD_150618.mat']);
NPDmat = fyA;
load([R.path.datapath '\nsPow_paper_RatNPD_150618.mat']);
nsPowmat = fyA;
load([R.path.datapath '\frq_paper_RatNPD_150618.mat']);
F_data = fxA(:,1);

meannpd_data = [];
% condsel = [1 2];
R.condnames = {'ON'};
R.Bcond = 1;
condsel = 1;
chsel = [1 3 2 4]'; % M1 STR GPe STN
klist = [2 3 1];

for C =1:numel(R.condnames)
    %     X = squeeze(fyA(:,:,:,:,2,:)); {i,j,dirc,cond,sub}
    for i = 1:size(chsel,1)
        for j = 1:size(chsel,1)
            if i==j
                F_model = R.frqz;
                % Log transform of group average
                Pxy = abs(log10(mean(vertcat(nsPowmat{chsel(i),condsel(C),:}),1)));
                % Put in the same frequency space as the models
                Pxy = interp1(F_data,Pxy,F_model');
                % Remove the powerline frequency
                Pxy(F_model>48 & F_model<52) = NaN;
                F_model(F_model>48 & F_model<52) = NaN;
                
                % Regress out the 1/f background
                [xCalc yCalc b Rsq] = linregress(log10(F_model)',Pxy);
                
                if nulldat == 1
                    Pxy = 10.^yCalc;
                else
                    Pxy = Pxy-yCalc;
                    % Bring to zero-base
                    Pxy = Pxy-min(Pxy);
                    Pxy = smoothdata(Pxy,'gaussian',round(4/diff(F_model(1:2))));
                    % Bring back to non-log space
                    Pxy = 10.^Pxy;
                end
                % Normalize and zero-base
                Pxy = (Pxy-mean(Pxy))./std(Pxy);
                Pxy = Pxy-min(Pxy);
%                 Pxy(F_model<12) = 0;
                meannpd_data(C,i,j,1,:) = Pxy;
            else
                
                for k = 1:size(NPDmat,3)
                    F_data = fxA(:,1);
                    F_model = R.frqz;
                    Cxy = mean(horzcat(NPDmat{chsel(i),chsel(j),klist(k),condsel(C),:})',1);
                    Cxy = Cxy.*1; % rescale 
                    
                    % Remove the powerline frequency
                    Cxy(F_model>48 & F_model<52) = [];
                    F_model(F_model>48 & F_model<52) = [];
                    % Put in the same frequency space as the models
                    Cxy = interp1(F_data,Cxy,F_model);
                    if nulldat == 1
                        Cxy = zeros(size(Cxy));
                    else
                    Cxy = smoothdata(Cxy','gaussian',round(4/diff(F_model(1:2))));
                    end
                    meannpd_data(C,i,j,k,:) = Cxy;
                end
                
            end
        end
    end
end

% Set data as working version
R.data.feat_emp{1} = meannpd_data;
% squeeze(meannpd_data(1,1,1,1,:))
R.data.feat_xscale{1} = R.frqz;

% Plot CSD
if plotop ==1
     R.plot.outFeatFx({R.data.feat_emp},[],R.data.feat_xscale,R,1,[])
end
