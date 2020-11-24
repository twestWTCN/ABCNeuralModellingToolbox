function [] = simAnPlotter_100717(feat_emp,feat_sim,F,R,bestn)
if isempty(feat_emp)
    feat_emp = {zeros(size(feat_sim{1}))};
end
if isempty(feat_sim)
    feat_sim= {zeros(size(feat_emp{1}))};
end
if nargin<4
    bestn = 1;
end
feat_data_n = feat_emp{1};
for L = 1:length(feat_sim)
    feat_sim_n = feat_sim{L};
    
    if L == bestn
        lwid = 2;
    else
        lwid = 0.5;
    end
    k = 0;
    
    if numel(size(feat_data_n)) > 2 % CSD data
        %% CSD PLOTS
        N = size(feat_data_n,1); M = size(feat_data_n,2);
        
        for i = 1:N
            for j = 1:M
                k = k+1;
                subplot(N,M,k)
                plot(F,squeeze(real(feat_sim_n(i,j,:))),'r--','linewidth',lwid); hold on
                plot(F,squeeze(real(feat_data_n(i,j,:))),'r','linewidth',2)
                
                plot(F,squeeze(imag(feat_sim_n(i,j,:))),'b--','linewidth',lwid)
                plot(F,squeeze(imag(feat_data_n(i,j,:))),'b','linewidth',2)
                
                plot(F,squeeze(abs(feat_sim_n(i,j,:))),'k--','linewidth',lwid)
                plot(F,squeeze(abs(feat_data_n(i,j,:))),'k','linewidth',2)
                
                xlim([min(R.frqz) max(R.frqz)])
                %         ylim([-0.03 0.03])
            end
        end
        set(gcf,'Position',[680         205        1112         893])
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
                    'LineStyle','none');
            end
        end
        %% TIME COURSE PLOTS
    else % time course
        N = size(feat_data_n,1);
        for i = 1:N
            subplot(1,N,i)
            plot(F,feat_sim_n(i,:),'r--','linewidth',lwid); hold on
            plot(F,feat_data_n(i,:),'r','linewidth',2)
            title(R.chloc_name{i})
             xlabel('Time'); ylabel(R.chloc_name{i})
        end
       
        set(gcf,'Position',[680         205        1112         893])
        
        
    end
end



