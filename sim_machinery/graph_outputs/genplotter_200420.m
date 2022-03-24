function genplotter_200420(datEmp,datSim,F,R,bestn,labelna)
if isempty(datEmp)
    datEmp = {zeros(size(datSim{1}))};
end
if ~isfield(R.plot,'featcolor')
    featcolor = [1 0 0];
else
    featcolor = R.plot.featcolor;
end
if isempty(datSim)
    for FN = 1:numel(R.data.datatype)
        datSim{1}{FN} = zeros(size(datEmp{1}{FN}));
    end
end
if nargin<5
    bestn = 1;
end
if isempty(labelna)
    labelna = 'NPD';
end
if ~isfield(R.plot,'holdop')
    R.plot.holdop = 0;
    clf
end
% Main Function Starts Here
for C = 1:numel(R.condnames)
    for FN = 1:numel(R.data.datatype)
        
        % hold on
        switch R.data.datatype{FN}
            case {'CSD','NPD'}
                figure(C*10 + FN)
                if ~R.plot.holdop; clf; end
                NPD_data_n = datEmp{1}{FN};
                
                for L = 1:length(datSim)
                    NPD_sim_n = datSim{L}{FN};
                    
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
                                plot(F{FN},squeeze(abs(NPD_sim_n(C,i,j,1,:))),'color',featcolor,'linestyle','--','linewidth',lwid); hold on
                                plot(F{FN},squeeze(imag(NPD_sim_n(C,i,j,1,:))),'b','linestyle','--','linewidth',lwid);
                            end
                            try
                                plot(F{FN},squeeze(abs(NPD_data_n(C,i,j,1,:))),'color',featcolor,'linewidth',2);
                                plot(F{FN},squeeze(imag(NPD_data_n(C,i,j,1,:))),'b','linewidth',2);
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
            case {'FANO','DURPDF','INTPDF'}
                figure(10 + FN)
                subplot(1,numel(R.condnames),C)
                if ~R.plot.holdop; cla; end
                
                fano_data = datEmp{1}{FN};
                plot(F{FN}(1:end),squeeze(fano_data(:,1,C)),'color',featcolor,'linewidth',2); hold on
                
                for L = 1:length(datSim)
                    fano_sim = datSim{L}{FN};
                    
                    if L == bestn
                        lwid = 2;
                    else
                        lwid = 0.5;
                    end
                    
                    plot(F{FN}(1:end),fano_sim(:,R.datinds,C),'color',featcolor,'linestyle','--','linewidth',lwid);
                end
                
                if strcmp(R.data.datatype{FN},'DURPDF')
                    xlabel('Burst duration (ms)');
                    ylabel('p.d.f')
                end
                
            case {'BRSTPROF','ENVPDF'}
                figure(10 + FN)
                subplot(1,numel(R.condnames),C)
                if ~R.plot.holdop; cla; end
                plot(F{FN},squeeze(datEmp{1}{FN}(:,:,C)),'color',featcolor,'linewidth',2); hold on
                for L = 1:length(datSim)
                    fano_sim = squeeze(datSim{L}{FN}(:,:,C));
                    
                    if L == bestn
                        lwid = 2;
                    else
                        lwid = 0.5;
                    end
                    
                    plot(F{FN},fano_sim(:,R.datinds),'color',featcolor,'linestyle','--','linewidth',lwid);
                end
                ylim([0 inf])
                if strcmp(R.data.datatype{FN},'BRSTPROF')
                    xlabel('Threshold');
                    ylabel('Mean Duration');
                elseif strcmp(R.data.datatype{FN},'ENVPDF')
                    xlabel('Burst amplitude');
                    ylabel('p.d.f')
                end
            otherwise
                warning('Plotting for datafeature not defined!')
                
        end
    end
end
