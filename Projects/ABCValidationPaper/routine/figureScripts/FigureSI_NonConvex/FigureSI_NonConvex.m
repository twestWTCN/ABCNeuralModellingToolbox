clear ; close all; %closeMessageBoxes
%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE SI- Non-Convexity of Objective Function
%%%%%%%%%%%%%%%%%%%%%%%%
% IF FRESH!
%      delete([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\MultiStartListWML.mat'])

%This should link to your repo folder
repopath = 'C:\Users\timot\Documents\GitHub\ABCNeuralModellingToolbox';
% repopath = 'C:\Users\Tim West\Documents\GitHub\ABCNeuralModellingToolbox'
addpath(repopath)
%This should be your projectname
projname = 'ABCValidationPaper';
R = ABCAddPaths(repopath,projname);

% Close all msgboxes
closeMessageBoxes

%% Setup for Model
R.out.tag = 'figureSI_NonConvex'; % Task tag
R = ABCsetup_partII_FullModel(R);

intag = 'figure5_ModelComp';
getParModulationData(R,intag)