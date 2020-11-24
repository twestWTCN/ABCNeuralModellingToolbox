function [] = comparplotter_180117(NPD_data,NPD_sim,F,R,bestn,labelna,leglist)
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
cmap = linspecer(length(NPD_sim));
for L = 1:length(NPD_sim)
    NPD_sim_n = NPD_sim{L};
    
lwid = 1;
    k = 0;
    N = size(NPD_data_n,1); M = size(NPD_data_n,2);
    for i = 1:N
        for j = 1:M
            k = k+1;
            subplot(N,M,k)
            if i == j
                plot(F{L},squeeze(NPD_sim_n(i,j,1,:)),'linewidth',lwid,'Color',cmap(L,:)); hold on
                                xlabel('Hz'); ylabel('Power'); title(sprintf('Ch %1.f Pxx',i))
%                                                     ylim([0 0.3])

            else
                a(L) = plot(F{L},squeeze(NPD_sim_n(i,j,2,:)),'linewidth',lwid,'Color',cmap(L,:)); hold on
                
%                 a(2) = plot(F{L},squeeze(NPD_sim_n(i,j,3,:)),'b--','linewidth',lwid)
                xlabel('Hz'); ylabel(labelna); title(sprintf('Ch %1.f -> Ch %1.f',i,j));
                if L==length(NPD_sim)
                legend(leglist)
                                    ylim([0 1.1])

                end
%                 plot(F,squeeze(NPD_sim_n(i,j,1,:)),'k--','linewidth',lwid)
%                 plot(F,squeeze(NPD_data_n(i,j,1,:)),'k','linewidth',2)
            end
            xlim([0 max(R.frqz)])
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
