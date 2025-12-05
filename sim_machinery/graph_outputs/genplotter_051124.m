% General Plotter for Empirical and Simulated Data
% This function generates plots comparing empirical and simulated data for different data types.
% Author: Timothy West, 2024, Imperial College London
% License: MIT License

function genplotter_051124(datEmp, datSim, F, R, bestn, labelna)
    % Set default values for input arguments
    if isempty(datEmp)
        for FN = 1:numel(R.data.datatype)
            datEmp{1}{FN} = zeros(size(datSim{1}{FN}));
        end
    end
    if ~isfield(R,'chsim_name')
        R.chsim_name = {};
    end
    if ~isfield(R.plot, 'featcolor')
        featcolor = {[1 0 0],[0 1 0]};
    else
        featcolor = {R.plot.featcolor(1,:),R.plot.featcolor(2,:)};
    end
    if isempty(datSim)
        for FN = 1:numel(R.data.datatype)
            datSim{1}{FN} = zeros(size(datEmp{1}{FN}));
        end
    end
    if nargin < 5
        bestn = 1;
    end
    if isempty(labelna)
        labelna = 'NPD';
    end
    if ~isfield(R.plot, 'holdop')
        R.plot.holdop = 0;
        cla;
    end

    % Determine subplot dimensions to create as square an array as possible
    if  numel(R.data.datatype)>1
        numPlots = numel(R.data.datatype);
        numRows = ceil(sqrt(numPlots));
        numCols = ceil(numPlots / numRows);
        DFLAG = 1; % multiplot flag
    elseif numel(R.chsim_name)>1
        numPlots =1;
        N = numel(R.chsim_name).^2;
        numRows = ceil(sqrt(N));
        numCols = ceil(N / numRows);
        DFLAG = 0; % multiplot flag
    else
        numPlots =1;
        N = 1;
        numRows = 1;
        numCols = 1;
        DFLAG = 0; % multiplot flag
    end

    % Main function starts here
    for C = 1:numel(R.condnames)
        figure(C);
        if ~R.plot.holdop; cla; end
        for FN = 1:numPlots
            if DFLAG == 1
                subplot(numRows, numCols, FN);
            end
            switch R.data.datatype{FN}
                case {'CSD', 'NPD'}
                    NPD_data_n = datEmp{1}{FN};
                    for L = 1:length(datSim)
                        NPD_sim_n = datSim{L}{FN};
                        lwid = (L == bestn) * 1.5 + 0.5; % Line width adjustment
                        [N, M] = size(NPD_sim_n, 2:3);
                        k = 0;
                        for i = 1:N
                            for j = 1:M
                                k = k + 1;
                                subplot(numRows, numCols, sub2ind([numCols numRows],j,i));
                                try
                                    plot(F{FN}, squeeze(abs(NPD_sim_n(C, i, j, 1, :))), 'color', featcolor{1}(i,:), 'linestyle', '--', 'linewidth', lwid); hold on;
                                    plot(F{FN}, squeeze(imag(NPD_sim_n(C, i, j, 1, :))), 'b', 'linestyle', '--', 'linewidth', lwid);
                                end
                                try
                                    plot(F{FN}, squeeze(abs(NPD_data_n(C, i, j, 1, :))), 'color', featcolor{2}(i,:), 'linewidth', 2);
                                    plot(F{FN}, squeeze(imag(NPD_data_n(C, i, j, 1, :))), 'b', 'linewidth', 2);
                                end
                                xlabel(R.plot.feat(FN).axtit{2});
                                ylabel(R.plot.feat(FN).axtit{1});
                                xlim(R.plot.feat(FN).axlim(1:2));
                                axis square;
                                if i == 1
                                    title(R.chsim_name{j});
                                elseif j == 1
                                    ylabel(R.chsim_name{i});
                                end
                            end
                        end
                    end
                case {'FANO', 'DURPDF', 'INTPDF'}
                    if ~R.plot.holdop; cla; end
                    fano_data = datEmp{1}{FN};
                    plot(F{FN}, squeeze(fano_data(:, 1, C)), 'color', featcolor{1}, 'linewidth', 2); hold on;
                    for L = 1:length(datSim)
                        fano_sim = datSim{L}{FN};
                        lwid = (L == bestn) * 1.5 + 0.5;
                        plot(F{FN}, fano_sim(:, R.datinds, C), 'color', featcolor{2}, 'linestyle', '--', 'linewidth', lwid);
                    end
                    xlabel(R.plot.feat(FN).axtit{2});
                    ylabel(R.plot.feat(FN).axtit{1});
                    xlim(R.plot.feat(FN).axlim(1:2));
                case {'BRSTPROF', 'ENVPDF'}
                    if ~R.plot.holdop; cla; end
                    plot(F{FN}, squeeze(datEmp{1}{FN}(:, :, C)), 'color', featcolor{1}, 'linewidth', 2); hold on;
                    for L = 1:length(datSim)
                        fano_sim = squeeze(datSim{L}{FN}(:, :, C));
                        lwid = (L == bestn) * 1.5 + 0.5;
                        plot(F{FN}, fano_sim(:, R.datinds), 'color', featcolor{2}, 'linestyle', '--', 'linewidth', lwid);
                    end
                    ylim([0, inf]);
                    xlabel(R.plot.feat(FN).axtit{2});
                    ylabel(R.plot.feat(FN).axtit{1});
                    xlim(R.plot.feat(FN).axlim(1:2));
                case 'time'
                    if ~R.plot.holdop; cla; end
                    LP = plot(F{FN}, datEmp{1}{FN}(:,:,C), 'linewidth', 2); hold on;
                    for i = 1:numel(LP)
                        LP(i).Color = featcolor{1}(i, :);
                    end
                    for L = 1:length(datSim)
                        xtmp = datSim{L}{FN}(:,:,C);
                        lwid = (L == bestn) * 1.5 + 0.5;
                        LP = plot(F{FN}, xtmp, 'linestyle', '--', 'linewidth', lwid);
                        for i = 1:numel(LP)
                            LP(i).Color = featcolor{2}(i, :);
                        end
                    end
                    xlabel(R.plot.feat(FN).axtit{2});
                    ylabel(R.plot.feat(FN).axtit{1});
                    xlim(R.plot.feat(FN).axlim(1:2));
                    % ylim(R.plot.feat(FN).axlim(3:4));
                otherwise
                    warning('Plotting for data feature not defined!');
            end
        end
    end
end
