% COMPUTE MODEL RELATIVE PROBABILITIES AND PLOT RESULTS
clear; close all
closeMessageBoxes

R = ABCAddPaths('Rat_NPD','rat_InDirect_ModelComp');
R.projectn = 'Rat_NPD';
R.out.tag = 'rat_InDirect_ModelComp'; % Task tag
R = simannealsetup_InDirect_ModelComp(R)

% Get empirical data
R = prepareRatData_InDirect_Group_NPD(R); 

%% Do the model probability computations
R.comptype = 1;
modelCompMaster(R,1:12) %,[1:8 10:12]

%% Plot the modComp results
R.modcomp.modN = [1:12];
R.modcompplot.NPDsel = [2 8 10]; %[6 9 10];
R.plot.confint = 'yes';
cmap = linspecer(numel(R.modcomp.modN));
cmap = cmap(end:-1:1,:);
plotModComp_091118(R,cmap)
figure(2)
subplot(3,1,1); ylim([-2 1])
