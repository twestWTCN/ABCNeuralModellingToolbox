%STN GPE CONFUSION MATRIX EVALUATION
%%%%%%%%%%%%%%%%%%%%%%%%
% rat_NPD_STN_GPe_Validation_ConfusionMat
% simAnnealAddPaths()
close all; clear
R = ABCAddPaths('Rat_NPD','rat_STN_GPe');

% Close all msgboxes
closeMessageBoxes()
rng(5353252)

%% Set Routine Pars
R.projectn = 'Rat_NPD'; % Project Name
R.out.tag = 'STN_GPe_ModComp'; % Task tag
R = simannealsetup_NPD_STN_GPe(R);

%% Set Routine Pars
load([R.rootn 'outputs\' R.out.tag '\ConfData'],'RSimData')
confmatlist = allcomb(1:3,1:3)';
R.tmp.confmat = confmatlist;
R.comptype = 2;
R.data = RSimData.data;
%% Do the model probability computations
% Create Job
% parallel.defaultClusterProfile('local')
% c = parcluster();
% c.NumWorkers = 4;
% j = createJob(c);
% for i = 1:6
%     createTask(j, @modelCompMaster,1,{R,1:9});
% end
% submit(j);wait(j);
% results = fetchOutputs(j);
% delete(j);

% modelCompMaster(R,1:9,[])
% modelCompMaster(R,2,[1 3])

%% Plot the modComp results
R.modcomp.modN = 1:9;
R.modcompplot.NPDsel = [1:9];
R.plot.confint = 'yes';
cmap = linspecer(numel(R.modcomp.modN));
plotModValidationEval(R,cmap)
figure(2)
subplot(3,1,1); ylim([-2 1])
