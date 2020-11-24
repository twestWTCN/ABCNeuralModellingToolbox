%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 2- Example inversion of STN/GPe Subcircuit
%%%%%%%%%%%%%%%%%%%%%%%%

clear ; close all
%This should link to your repo folder
repopath = 'C:\Users\timot\Documents\GitHub\ABCNeuralModellingToolbox';
%This should be your projectname
projname = 'ABCValidationPaper';
R = ABCAddPaths(repopath,projname);

% Close all msgboxes
closeMessageBoxes
rng(6439735)

%% Set Routine Pars
R.out.tag = 'figureSI_sensorNoise_abias'; % Task tag
R = ABCsetup_partI_STNGPe(R);
R.SimAn.pOptList = {'.int{src}.T','.C','.A','.D'}; %,'.S','.int{src}.G','.int{src}.S','.D','.A',,'.int{src}.BG','.int{src}.S','.S','.D','.obs.LF'};  %,'.C','.obs.LF'}; % ,'.obs.mixing','.C','.D',

%% Prepare the data
R = prepareRatData_STN_GPe_NPD(R);

%% Prepare Model
M = 1;
Clist = logspace(log10(0.001),log10(10),10);
for T = 1:10
    %     switch T
    %         case 1
    %             R.obs.gainmeth = {'unitvar'};
    %         case 2
    %             R.obs.gainmeth = {'obsnoise','unitvar'};
    %             R.obs.Cnoise = [0.02 0.02]; % Sensor Noise SNR prior i.e. 1/x signal to noise ratio
    %         case 3
    %             R.obs.gainmeth = {'obsnoise','unitvar'};
    %             R.obs.Cnoise = [0.2 0.2]; % Sensor Noise SNR prior i.e. 1/x signal to noise ratio
    %         case 4
    %             R.obs.gainmeth = {'obsnoise','unitvar'};
    %             R.obs.Cnoise = [2 2]; % Sensor Noise SNR prior i.e. 1/x signal to noise ratio
    %     end
    
    R.obs.Cnoise = [-1 Clist(T)]; %repmat(Clist(T),1,2);
    f = msgbox(sprintf('Probabilities for Model %.0f',T));
    
    modelspec = eval(['@MS_rat_STN_GPe_ModComp_Model' num2str(M)]);
    [R,p,m] = modelspec(R);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R.out.dag = sprintf([R.out.tag '_M%.0f'],T); % 'All Cross'
    R = setSimTime(R,28);
    R.Bcond = 0;
    R.plot.flag = 0;
    R.plot.save = 0;
    R.SimAn.convIt.dEps = 1e-8;
    [p] = SimAn_ABC_250320(R,p,m);
    closeMessageBoxes
end

%% Compare Model
R.comptype = 1;
modelCompMaster_160620(R,1:10,[])
R.modcomp.modN = 1:10;
R.modcompplot.NPDsel = [1:2:10];
R.plot.confint = 'yes';
R.plot.cmplx = 0;
cmap = linspecer(numel(R.modcomp.modN));
plotModComp_310520(R,cmap)
