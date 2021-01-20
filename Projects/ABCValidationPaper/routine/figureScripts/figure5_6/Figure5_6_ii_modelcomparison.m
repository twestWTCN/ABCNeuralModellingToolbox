clear ; close all;% closeMessageBoxes
%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 5/6- (II) MODEL COMPARISON
%%%%%%%%%%%%%%%%%%%%%%%%

%This should link to your repo folder
repopath = 'C:\Users\timot\Documents\GitHub\ABCNeuralModellingToolbox';
addpath(repopath)
%This should be your projectname
projname = 'ABCValidationPaper';
R = ABCAddPaths(repopath,projname);

%% Set Routine Pars
R.out.tag = 'figure5_ModelComp'; % Task tag
R = ABCsetup_partII_FullModel(R);

% IF FRESH START
% delete([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\WorkingModList.mat'])
%% Prepare the data
R = prepareRatData_InDirect_Group_NPD(R);

%% Do the model probability computations
R.comptype = 1;
R.out.tag2 = 'figure5_ModelComp';
% modelCompMaster_160620(R,1:12,[]) %,[1:8 10:12]
%% Plot the modComp results
R.modcomp.modN = [1:12];
R.modcompplot.NPDsel = [];% [10 8 4]; %[10 8 4]; %[8 10]; %[6 9 10];
R.plot.confint = 'yes';
cmap = linspecer(numel(R.modcomp.modN));
cmap = cmap(end:-1:1,:);
plotModCompABCValidationPaper(R,cmap)
% figure(2)
% subplot(3,1,1); ylim([-2 1])
