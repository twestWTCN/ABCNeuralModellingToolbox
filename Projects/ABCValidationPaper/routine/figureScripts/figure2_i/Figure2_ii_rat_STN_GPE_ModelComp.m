function Figure2_ii_rat_STN_GPE_ModelComp(R)
%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 2- Example model selection of STN/GPe Subcircuit
%%%%%%%%%%%%%%%%%%%%%%%%
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

