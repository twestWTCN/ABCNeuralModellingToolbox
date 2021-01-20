clear ; close all; closeMessageBoxes
%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 4- (I) Face Validation - Fitting Confusion Matrices with 8 node network
%%%%%%%%%%%%%%%%%%%%%%%%

%This should link to your repo folder
repopath = 'C:\Users\Tim West\Documents\GitHub\ABCNeuralModellingToolbox';
addpath(repopath)
%This should be your projectname
projname = 'ABCValidationPaper';
R = ABCAddPaths(repopath,projname);

R.out.tag = 'figure4_confusionMatrix';
R = ABCsetup_partII_FullModel(R);

%% First Simulate the Data from the Empirically Fitted Models
% delete([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\ConfWorkList.mat'])
% modlist = [11 3 8];
% for modID = 1:3
%     Rt = []; % Temp R Struc
%     % Recover Fitted Parameters
%     tagname = 'figure5_ModelComp';
%     [r2,pMAP{modlist(modID)},feat_sim{modlist(modID)},~,~,~,Rmod{modlist(modID)}] = getModelData(R,tagname,modlist(modID));
% end
% mkdir([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag])
% save([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\ConfData'],'feat_sim','pMAP','modlist','Rmod')
load([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\ConfData'],'feat_sim','pMAP','modlist','Rmod')
% Make matrix of combinations for confusion matrix
confmatlist = allcomb(modlist,modlist)';
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
        RSim = Rmod{SimMod};
        % Adjust R
        fprintf('Fitting Model %.0f to data %.0f',SimMod,SimData)
        f = msgbox(sprintf('Fitting Model %.0f to data %.0f',SimMod,SimData));
        
        modelspec = eval(['@MS_rat_InDrt_ModCompRev2_Model' num2str(SimMod)]);
        [~,p,m,uc] = modelspec(RSim);
%         pause(1)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        RSim.data.feat_emp = feat_sim{SimData};
        RSim.data.feat_xscale{1} = RSim.frqz;
        RSim.out.dag = sprintf([R.out.tag '_DataM%.0f_ParM%.0f'],SimData,SimMod); % 'All Cross'
        RSim.Bcond = 0;
        RSim.plot.flag = 1;
        RSim = rmfield(RSim,{'Mfit'}); % Important!! Reverts to unfitted model spec
        RSim.SimAn.convIt.dEps = 0.015;
         RSim = setSimTime(RSim,32);
        SimAn_ABC_201120(RSim,p,m);
        closeMessageBoxes
    end
end
