function FigurePLOS_LesionNPDOnly(R)
%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 5/6- (I) MODEL FITTING
%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Routine Pars
% R.out.tag = 'figurePLOS_NPDOnly'; % Task tag: NPDOnly and Lesion
R.out.tag = 'figurePLOS_Lesion'; % Task tag: NPDOnly and Lesion

R = ABCsetup_FigurePLOS(R);

% IF FRESH START
% delete([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\WorkingModList.mat'])
%% Prepare the data
R = prepareRatData_InDirect_Group_NPD(R);
modID = 10;

R.obs.LF = [1 1 1 1 1 1].*10; % Fit visually and for normalised data
R.obs.Cnoise = [0.2 0.2 0.2 0.2 0.2 0.2];
R.nmsim_name = {'MMC','STR','GPE','STN','GPI','THAL'};
R.chsim_name = {'MMC','STR','GPE','STN','GPI','THAL'};

fprintf('Now Fitting Model %.0f',modID)
f = msgbox(sprintf('Fitting Model %.0f',modID));

%% Prepare Model
modelspec = eval(['@MS_rat_InDrt_ModCompPLOSCBRev2_Model' num2str(modID)]);
[R p m uc] = modelspec(R); % M! intrinsics shrunk"
R.out.dag = sprintf([R.out.tag '_M%.0f'],modID);

%% Run ABC Optimization
R = setSimTime(R,32);
R.Bcond = 0;
R.plot.flag = 1;
% R.objfx.specspec = 'npd_only'; %
 [R,m,p] = loadABCData_160620(R);
[p] = SimAn_ABC_201120(R,p,m);
closeMessageBoxes

