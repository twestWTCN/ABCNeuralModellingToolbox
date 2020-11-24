%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 3- Multistarts
%%%%%%%%%%%%%%%%%%%%%%%%
% IF FRESH!
%  delete([R.rootn 'outputs\' R.out.tag '\MultiStartList.mat'])

clear ; close all
 closeMessageBoxes

%This should link to your repo folder
repopath = 'C:\Users\timot\Documents\GitHub\ABCNeuralModellingToolbox';
%This should be your projectname
projname = 'ABCValidationPaper';
R = ABCAddPaths(repopath,projname);

% Close all msgboxes
closeMessageBoxes
rng(6439735)

%% Set Routine Pars
R.out.tag = 'figure3_MultiStart'; % Task tag
R = ABCsetup_partI_STNGPe(R);


try
    load([R.path.rootn 'outputs\' R.out.tag '\MultiStartList'])
    disp('Loaded Mod List!!')
catch
    WML = [];
    mkdir([R.path.rootn 'outputs\' R.out.tag ]);
    save([R.path.rootn 'outputs\' R.out.tag '\MultiStartList'],'WML')
    disp('Making Mod List!!')
end
N = 10; % number of starts
%% Prepare Model
for multiStart = 1:2*N
    load([R.path.rootn 'outputs\' R.out.tag '\MultiStartList'],'WML')
    if ~any(intersect(WML,multiStart))
        WML = [WML multiStart];
        save([R.path.rootn 'outputs\' R.out.tag '\MultiStartList'],'WML')
        disp('Writing to Mod List!!')
        fprintf('Now Fitting Mulitstart %.0f',multiStart)
        f = msgbox(sprintf('Fitting Multistart %.0f',multiStart));
        
        if multiStart < N
            [~,pMAP{multiStart},feat_sim] = getModelData(R,'figure2_FitDemo',1);
            R.data.feat_emp = feat_sim;
            R.data.feat_xscale{1} = R.frqz;
        else
            [~,pMAP{multiStart},feat_sim] = getModelData(R,'figure2_FitDemo',2);
            R.data.feat_emp = feat_sim;
            R.data.feat_xscale{1} = R.frqz;
        end
        modelspec = eval(['@MS_rat_STN_GPe_ModComp_Model' num2str(1)]);
        [R,p,m] = modelspec(R);
        pause(2)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        R.out.dag = sprintf('NPD_STN_GPe_MultiStart_M%.0f',multiStart); % 'All Cross'
        % Delete Previous Saves
        delete([R.path.rootn 'outputs\' R.out.tag '\' R.out.dag '\*'])
        R.SimAn.rep = 128;
        R = setSimTime(R,32);
        R.Bcond = 0;
        R.plot.save = 1;
        R.plot.flag = 1;
        R.SimAn.convIt.eqN = 3;
        R.SimAn.convIt.dEps = 2e-3;
        R.SimAn.jitter = 1;
        SimAn_ABC_201120(R,p,m);
        closeMessageBoxes
    end
end

if i>inf
    modelspec = MS_rat_STN_GPe_ModComp_Model1(R);
    modelspec = MS_rat_STN_GPe_ModComp_Model2(R);
    modelspec = MS_rat_STN_GPe_ModComp_Model3(R);
end
