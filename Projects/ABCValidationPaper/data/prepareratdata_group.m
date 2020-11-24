function [nl] = prepareratdata_group(root,projectn)
% clear all
R = buildheader_rat();
srclist = {'M1','STN','GPe','STR'};
for cond = 1:2
    datacomb = [];
    for sub  = 1:length(R.subnames{cond})
        load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
        ftdata = FTdata.EpochedData;
        for i=1:4
            chind = min(find(strncmp(srclist{i},ftdata.label,length(srclist{i}))));
            chlist{i} = ftdata.label{chind};
            chind = [];
        end
        cfg = [];
        cfg.channel = chlist;
        ftdata = ft_selectdata(cfg,ftdata);
        chlist = [];
        datastch = [ftdata.trial{:}];
        datacomb = [datacomb datastch];
    end
    datastack{cond} =datacomb;
end

grandcat = [datastack{:}];
for i = 1:4
    grandmean(i) = mean(grandcat(i,:),2);
    grandstd(i) = std(grandcat(i,:),1);
end


for cond = 1:2
    datacomb = [];
    for sub  = 1:length(R.subnames{cond})
        load([R.analysispath R.pipestamp '\data\processed\' R.subnames{cond}{sub} '_' R.condnames{cond} '_' R.pipestamp '.mat'])
        for i=1:4
            chind = min(find(strncmp(srclist{i},ftdata.label,length(srclist{i}))));
            chlist{i} = ftdata.label{chind};
            chind = [];
        end
        cfg = [];
        cfg.channel = chlist;
        ftdata = ft_selectdata(cfg,ftdata);
        chlist = [];
        
        x = [];
        for j = 1:4
            for i = 1:size(ftdata.trial,2)
            dat = ftdata.trial{i}(j,:);
            dat = (dat-grandmean(j))./ grandstd(j);
            x(:,i,j) = dat;
            end
        end
        x = reshape(x,[],4);
        x= x';
        datacomb = [datacomb x];
        clear x
    end
    datastack{cond} =datacomb;
end



ftdata_control.label = {'MTX','STN','GPe','STR'}; %ftdata.label;
ftdata_control.trial{1} = datastack{1};
ftdata_control.fsample = 250;
ftdata_control.time = {linspace(0,length(datastack{1})/250,length(datastack{1}))};
% ftdata_control.dimord = 'chan_time';

ftdata_lesion.label = {'MTX','STN','GPe','STR'}; %ftdata.label;
ftdata_lesion.trial{1} = datastack{2};
ftdata_lesion.fsample = 250;
ftdata_lesion.time = {linspace(0,length(datastack{2})/250,length(datastack{2}))};

if length(ftdata_control.trial{1})> length(ftdata_lesion.trial{1})
    ftdata_control.trial{1} = ftdata_control.trial{1}(:,1:length(ftdata_lesion.trial{1})); ftdata_control.time = ftdata_lesion.time;
else
    ftdata_lesion.trial{1} = ftdata_lesion.trial{1}(:,1:length(ftdata_control.trial{1}));ftdata_lesion.time = ftdata_control.time;
end

mkdir([root],'data\storage')
filepathn = [root  '\data\storage']
save([filepathn '\average_rat_control'],'ftdata_control')
save([filepathn '\average_rat_lesion'],'ftdata_lesion')