function genplotter_200420(NPD_data,NPD_sim,F,R,bestn,labelna)
if isempty(NPD_data)
    NPD_data = {zeros(size(NPD_sim{1}))};
end
if ~isfield(R.plot,'cmap')
    R.plot.cmap = [1 0 0];
end
if isempty(NPD_sim)
    for FN = 1:numel(R.data.datatype)
    NPD_sim{1}{FN} = zeros(size(NPD_data{1}{FN}));
    end
end
if nargin<5
    bestn = 1;
end
if isempty(labelna)
    labelna = 'NPD';
end


% Main Function Starts Here


for FN = 1:numel(R.data.datatype)
    switch R.data.datatype{FN}
        case {'CSD','NPD'}
            NPD_data_n = NPD_data{1}{FN};
            O = size(NPD_data_n,1);

            for C = 1:O
                figure(C*10 + FN)
                clf
                for L = 1:length(NPD_sim)
                    NPD_sim_n = NPD_sim{L}{FN};
                    
                    if L == bestn
                        lwid = 2;
                    else
                        lwid = 0.5;
                    end
                    
                    N = size(NPD_sim_n,2); M = size(NPD_sim_n,3);
                    
                    k = 0;
                    for i = 1:N
                        for j = 1:M
                            k = k+1;
                            subplot(N,M,k)
                            try
                                plot(F{FN},squeeze(abs(NPD_sim_n(C,i,j,1,:))),'r--','linewidth',lwid); hold on
%                                 plot(F{FN},squeeze(imag(NPD_sim_n(C,i,j,1,:))),'b--','linewidth',lwid);
                            end
                            try
                                plot(F{FN},squeeze(abs(NPD_data_n(C,i,j,1,:))),'r','linewidth',2);
%                                 plot(F{FN},squeeze(imag(NPD_data_n(C,i,j,1,:))),'b','linewidth',2);
                            end
                            xlabel('Hz'); ylabel('Power'); %title(sprintf('Ch %1.f Pxx',i))
                            xlim([min(F{FN}) max(F{FN})])
                            axis square
                            %         ylim([-0.03 0.03])
                            if i==1
                                title(R.chsim_name{j})
                            elseif j == 1
                                ylabel(R.chsim_name{i})
                            end
                        end
                    end
                end
                set(gcf,'Position',[380         235        1112         893])
            end
            
        case {'FANO','DUR'}
            fano_data = NPD_data{1}{FN};
            figure(1*10 + FN)
            clf
            plot(F{2}(2:end),fano_data,'r','linewidth',2); hold on

            for L = 1:length(NPD_sim)
                fano_sim = NPD_sim{L}{FN};
                
                if L == bestn
                    lwid = 2;
                else
                    lwid = 0.5;
                end
                
                 plot(F{FN}(2:end),fano_sim(:,R.datinds),'r--','linewidth',lwid);
            end            
    end
end