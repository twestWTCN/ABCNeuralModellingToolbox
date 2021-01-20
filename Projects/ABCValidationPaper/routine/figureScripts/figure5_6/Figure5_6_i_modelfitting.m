clear ; close all;% closeMessageBoxes
%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 5/6- (I) MODEL FITTING
%%%%%%%%%%%%%%%%%%%%%%%%

%This should link to your repo folder
repopath = 'C:\Users\timot\Documents\GitHub\ABCNeuralModellingToolbox';
addpath(repopath)
%This should be your projectname
projname = 'ABCValidationPaper';
R = ABCAddPaths(repopath,projname);

%% Set Routine Pars
R.out.tag = 'figure5_ModelComp'; % Task tag
R = ABCsetup_partII_FullModel(R);

% IF FRESH START
% delete([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\WorkingModList.mat'])
%% Prepare the data
R = prepareRatData_InDirect_Group_NPD(R);
WML = [];

try
    load([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\WorkingModList'])
    disp('Loaded Mod List!!')
catch
    WML = [];
    mkdir([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag]);
    save([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\WorkingModList'],'WML')
    disp('Making Mod List!!')
end

for modID =1:12
    if modID>= 7
        R.obs.LF = [1 1 1 1 1 1].*10; % Fit visually and for normalised data
        R.obs.Cnoise = [0.2 0.2 0.2 0.2 0.2 0.2];
        R.nmsim_name = {'MMC','STR','GPE','STN','GPI','THAL'};
        R.chsim_name = {'MMC','STR','GPE','STN','GPI','THAL'};
    else
        R.chsim_name = {'MMC','STR','GPE','STN'};
        R.obs.LF = [1 1 1 1]*10; % Fit visually and for normalised data
        R.obs.Cnoise = [0.2 0.2 0.2 0.2];
    end
    load([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\WorkingModList'],'WML')
    if ~any(intersect(WML,modID))
        WML = [WML modID];
        save([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\WorkingModList'],'WML')
        disp('Writing to Mod List!!')
        fprintf('Now Fitting Model %.0f',modID)
        f = msgbox(sprintf('Fitting Model %.0f',modID));
        
        %% Prepare Model
        modelspec = eval(['@MS_rat_InDrt_ModCompRev2_Model' num2str(modID)]);
        [R p m uc] = modelspec(R); % M! intrinsics shrunk"
        pause(5)
        R.out.dag = sprintf([R.out.tag '_M%.0f'],modID);
        
        %% Run ABC Optimization
        R = setSimTime(R,32);
        R.Bcond = 0;
        R.plot.flag = 1;
        [p] = SimAn_ABC_201120(R,p,m);
        closeMessageBoxes
    end
end
