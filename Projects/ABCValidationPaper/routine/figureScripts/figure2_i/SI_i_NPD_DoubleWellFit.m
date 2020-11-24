%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 2- Example inversion of STN/GPe Subcircuit
%%%%%%%%%%%%%%%%%%%%%%%%

clear ; close all
%This should link to your repo folder
repopath = 'C:\Users\timot\Documents\GitHub\ABCNeuralModellingToolbox';
addpath(repopath)
%This should be your projectname
projname = 'ABCValidationPaper';
R = ABCAddPaths(repopath,projname);

% Close all msgboxes
closeMessageBoxes
rng('default'); rng(6439735)

%% Set Routine Pars
R.out.tag = 'SI_DoubleWell'; % Task tag
R = ABCsetup_SIDoubleWell(R);

%% Prepare the data
R = setSimTime(R,50);
R = prepareDoubleWellData(R);

%% Prepare Model
M = 1;
modelspec = eval(['@MS_doubleWellSpec_' num2str(M)]);
[R,p,m] = modelspec(R);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R.out.dag = sprintf([R.out.tag '_M%.0f'],M); % 'All Cross'
R.Bcond = 0;
R.plot.flag = 1;
R.plot.save = 1;
R.SimAn.convIt.dEps = 1e-12;
% R.SimAn.rep = 64;
SimAn_ABC_201120(R,p,m);
