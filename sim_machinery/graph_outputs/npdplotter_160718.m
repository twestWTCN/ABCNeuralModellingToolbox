function [] = npdplotter_160718(NPD_data,NPD_sim,F,R,bestn,labelna)
if isempty(NPD_data)
    NPD_data = {zeros(size(NPD_sim{1}))};
end
if isempty(NPD_sim)
    NPD_sim= {zeros(size(NPD_data{1}))};
end
if nargin<5
    bestn = 1;
end
if isempty(labelna)
    labelna = 'NPD';
end
NPD_data_n = NPD_data{1};
for L = 1:length(NPD_sim)
    NPD_sim_n = NPD_sim{L};
    
    if L == bestn
        lwid = 2;
    else
        lwid = 0.5;
    end
    k = 0;
    N = size(NPD_data_n,2); M = size(NPD_data_n,3); O = size(NPD_data_n,1);
    for C = 1:O
        k = 0;
        for i = 1:N
            for j = 1:M
                k = k+1;
                subplot(N,M,k)
                if i == j
                    plot(F,squeeze(NPD_sim_n(C,i,j,1,:)),'k--','linewidth',lwid); hold on
                    plot(F,squeeze(NPD_data_n(C,i,j,1,:)),'k','linewidth',2)
                    xlabel('Hz'); ylabel('Power'); %title(sprintf('Ch %1.f Pxx',i))
                else
                    a(1) = plot(F,squeeze(NPD_sim_n(C,i,j,2,:)),'r--','linewidth',lwid); hold on
                    plot(F,squeeze(NPD_data_n(C,i,j,2,:)),'r','linewidth',2)
                    
                    a(2) = plot(F,squeeze(NPD_sim_n(C,i,j,3,:)),'b--','linewidth',lwid);
                    plot(F,squeeze(NPD_data_n(C,i,j,3,:)),'b','linewidth',2)
                    xlabel('Hz'); ylabel(labelna); %title(sprintf('Ch %1.f -> Ch %1.f NPD',i,j));
                    %                 legend(a,{'Forward','Reverse'})
%                     hold on
%                     plot(F,squeeze(NPD_sim_n(C,i,j,1,:)),'k--','linewidth',lwid)
%                     plot(F,squeeze(NPD_data_n(C,i,j,1,:)),'k','linewidth',2)
                end
                xlim([min(R.frqz) max(R.frqz)])
                %         ylim([-0.03 0.03])
            end
        end
    end
end
set(gcf,'Position',[380         235        1112         893])
%LABEL LOCATIONS
% % Top row
% labPos(:,1,1) = [0.182654676258993 0.938409854423292 0.0548561151079137 0.0313549832026876];
% labPos(:,2,1) = [0.395100719424461 0.936326987681971 0.0548561151079137 0.0313549832026876];
% labPos(:,3,1) = [0.603266187050362 0.93552071668533 0.0548561151079137 0.0313549832026876];
% % labPos(:,4,1) = [0.7954964028777 0.935363941769317 0.0548561151079137 0.0313549832026876];
%
% % First Column
% labPos(:,1,2) =[0.0369352517985632 0.830414333706607 0.0548561151079137 0.0313549832026876];
% labPos(:,2,2) =[0.039633093525182 0.610929451287794 0.0548561151079137 0.0313549832026876];
% labPos(:,3,2) =[0.0425467625899302 0.392721164613662 0.0548561151079137 0.0313549832026876];
% % labPos(:,4,2) =[0.0436618705035993 0.178992161254199 0.0548561151079137 0.0313549832026876];
%
% for i = 1:N
%     for j = 1:M
%         annotation(gcf,'textbox',...
%             labPos(:,i,j)',...
%             'String',R.chloc_name{i},...
%             'LineStyle','none');
%     end
% end
%
%
