%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE SI- Examination of the role of sensor noise approximations
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
rng(6439735)

%% Set Routine Pars
R.out.tag = 'figureSI_sensorNoise_abias'; % Task tag
R = ABCsetup_partI_STNGPe(R);
R.SimAn.pOptList = {'.int{src}.T','.C','.A','.D'}; %,'.S','.int{src}.G','.int{src}.S','.D','.A',,'.int{src}.BG','.int{src}.S','.S','.D','.obs.LF'};  %,'.C','.obs.LF'}; % ,'.obs.mixing','.C','.D',
%% Prepare the data at different sensor noise levels
M = 1;
Clist = logspace(log10(0.001),log10(10),10);

%% Create Datasets
% for T =1:10
%         [~,pMAP{T},feat_sim{T}] = getModelData_sensorNoise(R,'figure2_FitDemo',1,repmat(Clist(T),1,2));
% 
% end
% save([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\SensorNoiseDataFeatures'],'pMAP','feat_sim')
load([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\SensorNoiseDataFeatures'],'pMAP','feat_sim')





%% Prepare Model
for T = 2:10
se = [-1 Clist(T)]; %repmat(Clist(T),1,2);
    f = msgbox(sprintf('Probabilities for Model %.0f',T));
    modelspec = eval(['@MS_rat_STN_GPe_ModComp_Model' num2str(M)]);
    [R,p,m] = modelspec(R);
    R.data.feat_emp = feat_sim{T};
    R.data.feat_xscale{1} = R.frqz;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R.out.dag = sprintf([R.out.tag '_M%.0f'],T); % 'All Cross'
    R = setSimTime(R,32);
    R.Bcond = 0;
    R.plot.flag = 0;
    R.plot.save = 0;
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
