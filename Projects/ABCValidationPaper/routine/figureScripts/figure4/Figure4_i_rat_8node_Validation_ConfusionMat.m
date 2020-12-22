clear ; close all; closeMessageBoxes
%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 4- (I) Face Validation - Confusion Matrices with 8 node network
%%%%%%%%%%%%%%%%%%%%%%%%

%This should link to your repo folder
repopath = 'C:\Users\timot\Documents\GitHub\ABCNeuralModellingToolbox';
addpath(repopath)
%This should be your projectname
projname = 'ABCValidationPaper';
R = ABCAddPaths(repopath,projname);

R.out.tag = 'figure4_confusionMatrix';

%% First Simulate the Data from the Empirically Fitted Models
delete([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\ConfWorkList.mat'])
for modID = 1:3
    Rt = []; % Temp R Struc
    % Recover Fitted Parameters
    tagname = 'figure5_ModelComp';
    [r2,pMAP{modID},feat_sim{modID}] = getModelData(R,tagname,modID);
end
save([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\ConfData'],'RSimData')
load([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\ConfData'],'RSimData')
% Make matrix of combinations for confusion matrix
confmatlist = allcomb(1:3,1:3)';
R.tmp.confmat = confmatlist;

% Create List (for parallelization across multiple MATLAB instances)
try
    load([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\ConfWorkList'])
    disp('Loaded Mod List!!')
catch
    WML = [];
    mkdir([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag ]);
    save([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\ConfWorkList'],'WML')
    disp('Making Mod List!!')
end

%% Prepare Model
for i =1:size(confmatlist,2)
    load([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\ConfWorkList'],'WML')
    if ~any(intersect(WML,i))
        WML = [WML i];
        save([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\ConfWorkList'],'WML')
        disp('Writing to Mod List!!')
        
        SimData = confmatlist(1,i);
        SimMod = confmatlist(2,i);
        
        fprintf('Fitting Model %.0f to data %.0f',SimMod,SimData)
        f = msgbox(sprintf('Fitting Model %.0f to data %.0f',SimMod,SimData));
        
        modelspec = eval(['@MS_rat_STN_GPe_ModComp_Model' num2str(SimMod)]);
        [R p m uc] = modelspec(R);
        pause(5)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        R.data.feat_emp = feat_sim{2};
        R.data.feat_xscale{1} = R.frqz;
        R.out.dag = sprintf([R.out.tag '_DataM%.0f_ParM%.0f'],SimData,SimMod); % 'All Cross'
        R.Bcond = 0;
        R.plot.flag = 1;
        SimAn_ABC_201120(R,p,m);
        closeMessageBoxes
    end
end
