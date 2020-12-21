%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 2- Example model selection of STN/GPe Subcircuit
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
rng(6439735)

%% Set Routine Pars
R.out.tag = 'figure2_FitDemo'; % Task tag
R = ABCsetup_partI_STNGPe(R);

%% Prepare the data

%% Do the model probability computations
R.comptype = 1;
modelCompMaster_160620(R,1,[]);

%% Plot the modComp results
R.modcomp.modN = 1:3;
R.modcompplot.NPDsel = [1:3];
R.plot.confint = 'yes';
cmap = linspecer(numel(R.modcomp.modN));
computeModComp_231120(R,cmap)


plotModComp_091118(R,cmap)
figure(2)
subplot(4,1,1); ylim([-2 1])
% subplot(4,1,2); ylim([0 3])

