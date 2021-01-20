function  [hl, hp, dl,flag,CSD_mean] = PlotFeatureConfInt_gen170620(R,permMod,fighan,cmap)
if ~isfield(R.plot,'confint')
    R.plot.confint = 1; % If not specified then dont plot confidence intervals
end
if ~isfield(R.plot,'cmplx')
    R.plot.cmplx = 0;
end
for featN = 1:numel(R.data.datatype)
    switch R.data.datatype{featN}
        case {'CSD','NPD'}
            nanflag = 1; pl = 1;
            while nanflag == 1
                CSD_data_n = permMod.feat_rep{pl}{featN};
                if any(isnan(CSD_data_n(:)))
                    pl = pl + 1;
                else
                    nanflag = 0;
                end
            end
            list =find(~isinf(permMod.r2rep) & ~isnan(permMod.r2rep));
            C = size(CSD_data_n,1); N = size(CSD_data_n,2); M = size(CSD_data_n,3);O = size(CSD_data_n,4);
            for cond = 1:C
                for ii = 1:size(list,2)
                    for i = 1:N
                        for j = 1:M
                            for k = 1:O
                                CSD_bank(:,i,j,k,ii) = permMod.feat_rep{list(ii)}{featN}(cond,i,j,k,:);
                            end
                        end
                    end
                end
                
                if nargin<3
                    figure
                else
                    figure(fighan(featN,cond))
                end
                
                if strncmp(R.data.datatype{featN},'CSD',3)
                    partlabs ={'Abs','Imag','Real'}; msr = 'CSD';
                    
                    CSD_mean = prctile(CSD_bank,50,5);
                    CSD_std(:,:,:,:,1)  = prctile(CSD_bank,50,5)-prctile(CSD_bank,25,5); %std(CSD_bank,1,5); %prctile(CSD_bank,5,5); %
                    CSD_std(:,:,:,:,2)  = prctile(CSD_bank,75,5)-prctile(CSD_bank,50,5); %std(CSD_bank,1,5); %prctile(CSD_bank,5,5); %
                elseif strncmp(R.data.datatype{featN},'NPD',3)
                    partlabs ={'Instant','Forward','Backward'}; msr = 'NPD';
                    CSD_mean = prctile(CSD_bank,50,5);
                    CSD_std(:,:,:,:,1)  = prctile(CSD_bank,50,5)-prctile(CSD_bank,25,5); %std(CSD_bank,1,5); %prctile(CSD_bank,5,5); %
                    CSD_std(:,:,:,:,2)  = prctile(CSD_bank,75,5)-prctile(CSD_bank,50,5); %std(CSD_bank,1,5); %prctile(CSD_bank,5,5); %
                end
                % F = repmat(R.frqz',1,3);
                F = linspace(min(R.frqz),max(R.frqz),size(CSD_mean,1));
                F = repmat(F',1,3);
                k = 0;
                for i = 1:N
                    for j = 1:M
                        k = k+1;
                        if i==j
                            c = 1; lr = 1;
                            q = 1;
                        elseif i>j
                            c = 2; lr=2;
                            q = 2;
                        elseif j>i
                            c = 3; lr=2;
                            q = 2;
                        end
                        
                        subplot(N,M,k)
                        Y = squeeze(CSD_mean(:,i,j,:));
                        B = squeeze(CSD_std(:,i,j,:,:));
                        B = permute(B,[1 3 2]);
                        alpval = 0.45;
                        if isequal(R.plot.confint,'none')
                            B = zeros(size(B));
                            alpval = 0;
                        end
                        if R.plot.cmplx
                            [hl, hp] = boundedline(F(:,lr),real(Y(:,lr)),real(B(:,:,lr)),'cmap',cmap,'alpha','transparency',alpval);
                            hl(1).LineWidth = 2; %hl(2).LineWidth = 1; %hl(3).LineWidth = 1;
                            hl(1).LineStyle = '-';% hl(2).LineStyle = '--';% hl(3).LineStyle = '--';
                            
                            [hl, hp] = boundedline(F(:,lr),imag(Y(:,lr)),imag(B(:,:,lr)),'cmap',cmap,'alpha','transparency',alpval);
                            hl(1).LineWidth = 2; %hl(2).LineWidth = 1; %hl(3).LineWidth = 1;
                            hl(1).LineStyle = '--';% hl(2).LineStyle = '--';% hl(3).LineStyle = '--';
                            if ~isequal(R.plot.confint,'none')
                                hout = outlinebounds(hl, hp);
                                set(hout(1),'LineWidth',1,'LineStyle','--'); %set(hout(2),'LineWidth',1);% set(hout(3),'LineWidth',1);
                            end
                            
                            hold on
                            for L = lr
                                dl = plot(R.data.feat_xscale{featN},real(squeeze(R.data.feat_emp{featN}(cond,i,j,L,:))),'color',[0 0 0],'LineWidth',1.5); hold on
                                dl = plot(R.data.feat_xscale{featN},imag(squeeze(R.data.feat_emp{featN}(cond,i,j,L,:))),'color',[0 0 0],'LineWidth',1.5,'LineStyle','--'); hold on
                            end
                            
                            
                        else
                            [hl, hp] = boundedline(F(:,lr),abs(Y(:,lr)),abs(B(:,:,lr)),'cmap',cmap,'alpha','transparency',alpval);
                            hl(1).LineWidth = 2; %hl(2).LineWidth = 1; %hl(3).LineWidth = 1;
                            hl(1).LineStyle = '-';% hl(2).LineStyle = '--';% hl(3).LineStyle = '--';
                            if ~isequal(R.plot.confint,'none')
                                hout = outlinebounds(hl, hp);
                                set(hout(1),'LineWidth',1,'LineStyle','--'); %set(hout(2),'LineWidth',1);% set(hout(3),'LineWidth',1);
                            end
                            
                            hold on
                            for L = lr
                                dl = plot(R.data.feat_xscale{featN},q.*abs(squeeze(R.data.feat_emp{featN}(cond,i,j,L,:))),'color',[0 0 0],'LineWidth',1.5); hold on
                            end
                        end
                        if i == j
                            ylim([0 5])
                            ylabel('Power')
                            title(R.chsim_name{i})
                        else
                            %                 ylim([0 0.5])
                            ylabel(msr)
                        end
                        xlim([min(R.frqz) max(R.frqz)])
                        xlabel('Freq (Hz)');
                        grid on
                        %         ylim([-0.03 0.03])
                    end
                end
                flag = 0;
                set(gcf,'Position',[680         112        1112         893])
            end
        case {'FANO','DUR'}
            fanodata = permMod.feat_rep{1}{featN};
            list =find(~isinf(permMod.ACCrep) & ~isnan(permMod.ACCrep) & (permMod.ACCrep>prctile(permMod.ACCrep,75)));
            N = size(fanodata,2);
            fano_std = []; fanobank = [];
            for ii = 1:size(list,2)
                for i = 1:N
                    fanobank(:,i,ii) = permMod.feat_rep{list(ii)}{featN}(:,i);
                    
                end
            end
            figure(fighan(featN,1))
            
            fanomean = prctile(fanobank,50,3);
            fano_std(:,:,1)  = std(fanobank,1,3); %prctile(fanobank,50,3)-prctile(fanobank,25,3); % %prctile(CSD_bank,5,5); %
            fano_std(:,:,2)  = std(fanobank,1,3); %prctile(fanobank,75,3)-prctile(fanobank,50,3); %std(CSD_bank,1,5); %prctile(CSD_bank,5,5); %
            
            F =  R.data.feat_xscale{featN};
            
            [hl, hp] = boundedline(F(2:end),fanomean(:,1),squeeze(fano_std(:,1,:)),'cmap',cmap,'alpha','transparency',alpval);
            hl(1).LineWidth = 2; %hl(2).LineWidth = 1; %hl(3).LineWidth = 1;
            hl(1).LineStyle = '-';% hl(2).LineStyle = '--';% hl(3).LineStyle = '--';
            
            hold on
            dl = plot(F(2:end),R.data.feat_emp{featN} ,'color',[0 0 0],'LineWidth',1.5); hold on
            
            xlabel('Burst Duration (ms)'); ylabel('p.d.');
            grid on
    end
    
end
% % %LABEL LOCATIONS
% % % Top row
% % labPos(:,1,1) = [0.182654676258993 0.938409854423292 0.0548561151079137 0.0313549832026876];
% % labPos(:,2,2) = [0.395100719424461 0.936326987681971 0.0548561151079137 0.0313549832026876];
% % labPos(:,3,3) = [0.603266187050362 0.93552071668533 0.0548561151079137 0.0313549832026876];
% % labPos(:,4,4) = [0.7954964028777 0.935363941769317 0.0548561151079137 0.0313549832026876];
% %
% % % First Column
% % labPos(:,1,1) =[0.0369352517985632 0.830414333706607 0.0548561151079137 0.0313549832026876];
% % labPos(:,2,2) =[0.039633093525182 0.610929451287794 0.0548561151079137 0.0313549832026876];
% % labPos(:,3,3) =[0.0425467625899302 0.392721164613662 0.0548561151079137 0.0313549832026876];
% % labPos(:,4,4) =[0.0436618705035993 0.178992161254199 0.0548561151079137 0.0313549832026876];
% %
% % for i = 1:N
% %     for j = 1:M
% %         annotation(gcf,'textbox',...
% %             labPos(:,i,j)',...
% %             'String',R.chsim_name{i},...
% %             'FontWeight','Bold',...
% %             'LineStyle','none');
% %     end
% % end
% % annotation(gcf,'textbox',...
% %     [0.24 0.013 0.54 0.05],...
% %     'String',{'Mean, Upper and Lower Quartile of the Posterior Distribution for the',' ABC Optimized Model Fit to Lesion Group Average CSD'},...
% %     'HorizontalAlignment','center',...
% %     'FitBoxToText','off',...
% %     'FontWeight','Bold',...,
% %     'FontSize',12,...
% %     'LineStyle','none');
% % % h = legend({partlabs{1},partlabs{2},partlabs{3},[partlabs{1} ' Mean'],[partlabs{2} ' Mean'],[partlabs{3}  ' mean'],[partlabs{1} ' Qrt.'],[partlabs{2} ' Qrt.'],[partlabs{3} ' Qrt']});
% % % set(h,...
% % %     'Position',[0.87 0.013 0.11 0.17]);
