function PlotFeatureConfInt_gen(R,d)
if nargin<2
    d = sprintf('%d',[R.d(1:3)]);
end

% load([R.rootn 'outputs\' R.out.tag '2\permMod_' R.out.tag '_' d '.mat'],'permMod')
load('permod.mat')
CSD_data_n = permMod.feat_rep{1};
N = size(CSD_data_n,1); M = size(CSD_data_n,2);
for ii = 1:size(permMod.feat_rep,2)
    for i = 1:N
        for j = 1:M
            for k = 1:O
                CSD_bank(:,i,j,ii) = permMod.feat_rep{ii}(i,j,:);
            end
        end
    end
end
CSD_mean(:,:,:,1) = mean(real(CSD_bank),4);
% CSD_std(:,:,:,1)  = std(real(CSD_bank),1,4);]
CSD_std(:,:,:,1,1)  = prctile(real(CSD_bank),25,4);
CSD_std(:,:,:,2,1)  = prctile(real(CSD_bank),75,4);

CSD_mean(:,:,:,2) = mean(imag(CSD_bank),4);
% CSD_std(:,:,:,2)  = std(imag(CSD_bank),1,4);
CSD_std(:,:,:,1,2)  = prctile(imag(CSD_bank),25,4);
CSD_std(:,:,:,2,2)  = prctile(imag(CSD_bank),75,4);

CSD_mean(:,:,:,3) = mean(abs(CSD_bank),4);
% CSD_std(:,:,:,3)  = std(abs(CSD_bank),1,4);
CSD_std(:,:,:,1,3)  = prctile(abs(CSD_bank),25,4);
CSD_std(:,:,:,2,3)  = prctile(abs(CSD_bank),75,4);

% F = repmat(R.frqz',1,3);
F = linspace(min(R.frqz),max(R.frqz),size(CSD_mean,1));
F = repmat(F',1,3);
cmap = linspecer(2);
cmap = [0 0 0; cmap];
 k = 0;
for i = 1:N
    for j = 1:M
        k = k+1;
        subplot(N,M,k)
        Y = squeeze(CSD_mean(:,i,j,:));
        B = squeeze(CSD_std(:,i,j,:,:));
        [hl, hp] = boundedline(F,Y,B,'cmap',cmap,'alpha','transparency',0.45);
        hl(1).LineWidth = 2; hl(2).LineWidth = 1; hl(3).LineWidth = 1;
        hout = outlinebounds(hl, hp);
        set(hout(1:3),'LineWidth',1,'LineStyle','--'); set(hout(2),'LineWidth',1); set(hout(3),'LineWidth',1);
        if i == j
            ylim([0 0.06])
        else
%             ylim([-6e-3 35e-3])
        end
        xlim([min(R.frqz) max(R.frqz)])
        xlabel('Freq (Hz'); ylabel('Amplitude')
        %         ylim([-0.03 0.03])
    end
end
    
    set(gcf,'Position',[680         112        1112         893])
%LABEL LOCATIONS
% Top row
labPos(:,1,1) = [0.182654676258993 0.938409854423292 0.0548561151079137 0.0313549832026876];
labPos(:,2,1) = [0.395100719424461 0.936326987681971 0.0548561151079137 0.0313549832026876];
labPos(:,3,1) = [0.603266187050362 0.93552071668533 0.0548561151079137 0.0313549832026876];
labPos(:,4,1) = [0.7954964028777 0.935363941769317 0.0548561151079137 0.0313549832026876];

% First Column
labPos(:,1,2) =[0.0369352517985632 0.830414333706607 0.0548561151079137 0.0313549832026876];
labPos(:,2,2) =[0.039633093525182 0.610929451287794 0.0548561151079137 0.0313549832026876];
labPos(:,3,2) =[0.0425467625899302 0.392721164613662 0.0548561151079137 0.0313549832026876];
labPos(:,4,2) =[0.0436618705035993 0.178992161254199 0.0548561151079137 0.0313549832026876];

for i = 1:4
    for j = 1:2
        annotation(gcf,'textbox',...
            labPos(:,i,j)',...
            'String',R.chloc_name{i},...
            'FontWeight','Bold',...
            'LineStyle','none');
    end
end

annotation(gcf,'textbox',...
    [0.24 0.013 0.54 0.05],...
    'String',{'Mean, Upper and Lower Quartile of the Posterior Distribution for the',' ABC Optimized Model Fit to Lesion Group Average CSD'},...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'FontWeight','Bold',...,
    'FontSize',12,...
    'LineStyle','none');
h = legend({'Abs','Imag','Real','Abs Mean','Real Mean','Real mean','Abs Qrt.','Imag Qrt.','Real Qrt'});
set(h,...
    'Position',[0.87 0.013 0.11 0.17]);
