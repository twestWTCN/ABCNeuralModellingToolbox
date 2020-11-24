function [] = csdplotter_220517(CSD_data,CSD_sim,F,R,bestn)
if isempty(CSD_data)
    CSD_data = {zeros(size(CSD_sim{1}))};
end
if isempty(CSD_sim)
    CSD_sim= {zeros(size(CSD_data{1}))};
end
if nargin<6
    bestn = 1;
end
CSD_data_n = CSD_data{1};
for L = 1:length(CSD_sim)
    CSD_sim_n = CSD_sim{L};
    
    if L == bestn
        lwid = 2;
    else
        lwid = 0.5;
    end
    k = 0;
    N = size(CSD_data_n,1); M = size(CSD_data_n,2);
    for i = 1:N
        for j = 1:M
            k = k+1;
            subplot(N,M,k)
            hold on
            if i==j
                plot(F,squeeze(abs(CSD_sim_n(i,j,:))),'k--','linewidth',lwid)
                plot(F,squeeze(abs(CSD_data_n(i,j,:))),'k','linewidth',2)
            else
                plot(F,squeeze(real(CSD_sim_n(i,j,:))),'r--','linewidth',lwid); hold on
                plot(F,squeeze(real(CSD_data_n(i,j,:))),'r','linewidth',2)
                
                plot(F,squeeze(imag(CSD_sim_n(i,j,:))),'b--','linewidth',lwid)
                plot(F,squeeze(imag(CSD_data_n(i,j,:))),'b','linewidth',2)
            end
            
            xlim([min(R.frqz) 45]) %max(R.frqz)])
            %         ylim([-0.03 0.03])
        end
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


