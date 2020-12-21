%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 2- Example inversion of STN/GPe Subcircuit
%%%%%%%%%%%%%%%%%%%%%%%%

clear ; close all
%This should link to your repo folder
% repopath = 'C:\Users\timot\Documents\GitHub\ABCNeuralModellingToolbox';
repopath = 'C:\Users\Tim West\Documents\GitHub\ABCNeuralModellingToolbox';

addpath(repopath)
%This should be your projectname
projname = 'ABCValidationPaper';
R = ABCAddPaths(repopath,projname);

% Close all msgboxes
closeMessageBoxes
rng('default'); rng(6439735)

%% Set Routine Pars
R.out.tag = 'figure2_FitDemo'; % Task tag
R = ABCsetup_partI_STNGPe(R);

%% Prepare the data
R = prepareRatData_STN_GPe_NPD(R);

%% Prepare Model
for M = 1:3
modelspec = eval(['@MS_rat_STN_GPe_ModComp_Model' num2str(M)]);
[R,p,m] = modelspec(R);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R.out.dag = sprintf([R.out.tag '_M%.0f'],M); % 'All Cross'
R = setSimTime(R,28);
R.Bcond = 0;
R.plot.flag = 1;
R.plot.save = 0;
[p] = SimAn_ABC_201120(R,p,m);
end
