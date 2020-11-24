function PlotFeatureConfInt_gen(R)
d =  R.out.dag;
% if nargin<2
%     d = sprintf('%d',[R.d(1:3)]);
% end
% clear all; close all
load([R.rootn 'outputs\' R.out.tag '2\permMod_' R.out.tag '_' d '.mat'],'permMod')

CSD_data_n = permMod.feat_rep{1};
% list = find([permMod.r2rep{:}]>prctile([permMod.r2rep{:}],1));
list = find([permMod.r2rep{:}]>prctile([permMod.r2rep{:}],85));

C = size(CSD_data_n,1); N = size(CSD_data_n,2); M = size(CSD_data_n,3);O = size(CSD_data_n,4);

for ii = 1:size(list,2)
    for i = 1:N
        for j = 1:M
            for k = 1:O
                CSD_bank(:,i,j,k,ii) = permMod.feat_rep{list(ii)}(C,i,j,k,:);
            end
        end
    end
end
if strncmp(R.data.datatype,'CSD',3)
    partlabs ={'Abs','Imag','Real'}; msr = 'CSD';
    
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
elseif strncmp(R.data.datatype,'NPD',3)
    partlabs ={'Instant','Forward','Backward'}; msr = 'NPD';
    CSD_mean = mean(CSD_bank,5);
    CSD_std(:,:,:,:,1)  = std(CSD_bank,1,5);%./sqrt(size(CSD_bank,5)); %prctile(CSD_bank,25,5);
    CSD_std(:,:,:,:,2)  = std(CSD_bank,1,5);%./sqrt(size(CSD_bank,5)); %prctile(CSD_bank,75,5);
end
% F = repmat(R.frqz',1,3);
F = linspace(min(R.frqz),max(R.frqz),size(CSD_mean,1));
F = repmat(F',1,3);
cmap = linspecer(2);
figure
cmap = [0 0 0; cmap];
k = 0;
for i = 1:N
    for j = 1:M
        k = k+1;
        if i==j
            c = 1; lr = 1;
        elseif i>j
            c = 2; lr=2;
        elseif j>i
            c = 3; lr=2;
        end        
        
        subplot(N,M,k)
        Y = squeeze(CSD_mean(:,i,j,:));
        B = squeeze(CSD_std(:,i,j,:,:));
        B = permute(B,[1 3 2]);
%         if i==j
%             B = B./10; %size(CSD_bank,5);
%         end
        [hl, hp] = boundedline(F(:,lr),Y(:,lr),B(:,lr),'cmap',cmap([c],:),'alpha','transparency',0.45);
        hl(1).LineWidth = 2; %hl(2).LineWidth = 1; %hl(3).LineWidth = 1; 
        hl(1).LineStyle = '--';% hl(2).LineStyle = '--';% hl(3).LineStyle = '--';
        hout = outlinebounds(hl, hp);
        set(hout(1),'LineWidth',1,'LineStyle','--'); %set(hout(2),'LineWidth',1);% set(hout(3),'LineWidth',1);
        hold on
        for L = lr
         plot(R.data.feat_xscale,squeeze(R.data.feat_emp(C,i,j,L,:)),'color',cmap(c,:)); hold on
        end
        if i == j
            ylim([0 4])
            ylabel('Amplitude')
        else
            ylim([0 0.3])
            ylabel(msr)
        end
        xlim([min(R.frqz) max(R.frqz)])
        xlabel('Freq (Hz)'); 
        %         ylim([-0.03 0.03])
    end
end

set(gcf,'Position',[680         112        1112         893])
%LABEL LOCATIONS
% Top row
labPos(:,1,1) = [0.182654676258993 0.938409854423292 0.0548561151079137 0.0313549832026876];
labPos(:,2,2) = [0.395100719424461 0.936326987681971 0.0548561151079137 0.0313549832026876];
labPos(:,3,3) = [0.603266187050362 0.93552071668533 0.0548561151079137 0.0313549832026876];
labPos(:,4,4) = [0.7954964028777 0.935363941769317 0.0548561151079137 0.0313549832026876];

% First Column
labPos(:,1,1) =[0.0369352517985632 0.830414333706607 0.0548561151079137 0.0313549832026876];
labPos(:,2,2) =[0.039633093525182 0.610929451287794 0.0548561151079137 0.0313549832026876];
labPos(:,3,3) =[0.0425467625899302 0.392721164613662 0.0548561151079137 0.0313549832026876];
labPos(:,4,4) =[0.0436618705035993 0.178992161254199 0.0548561151079137 0.0313549832026876];

for i = 1:N
    for j = 1:M
        annotation(gcf,'textbox',...
            labPos(:,i,j)',...
            'String',R.chsim_name{i},...
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
% h = legend({partlabs{1},partlabs{2},partlabs{3},[partlabs{1} ' Mean'],[partlabs{2} ' Mean'],[partlabs{3}  ' mean'],[partlabs{1} ' Qrt.'],[partlabs{2} ' Qrt.'],[partlabs{3} ' Qrt']});
% set(h,...
%     'Position',[0.87 0.013 0.11 0.17]);
