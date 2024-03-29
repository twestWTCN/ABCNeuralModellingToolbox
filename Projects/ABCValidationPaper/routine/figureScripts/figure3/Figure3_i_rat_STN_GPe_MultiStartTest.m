function Figure3_i_rat_STN_GPe_MultiStartTest(R)
%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 3- (I) Multistart Fitting
%%%%%%%%%%%%%%%%%%%%%%%%
% Close all msgboxes
closeMessageBoxes

%% Set Routine Pars
R.out.tag = 'figure3_MultiStart'; % Task tag
R = ABCsetup_partI_STNGPe(R);
R.SimAn.pOptList = {'.int{src}.T','.int{src}.S','.C','.A','.D'}; %,'.S','.int{src}.G','.int{src}.S','.D','.A',,'.int{src}.BG','.int{src}.S','.S','.D','.obs.LF'};  %,'.C','.obs.LF'}; % ,'.obs.mixing','.C','.D',

%% Create Datasets
for DS = 1:2
    if DS == 1
        [~,pMAP{DS},feat_sim{DS}] = getMultiStartData(R,'figure2_FitDemo',1);
    else
        [~,pMAP{DS},feat_sim{DS}] = getMultiStartData(R,'figure2_FitDemo',2);
    end
end
save([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\MultiStartDataFeatures'],'pMAP','feat_sim')
% load([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\MultiStartDataFeatures'],'pMAP','feat_sim')

try
    load([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\MultiStartListWML'])
    disp('Loaded Mod List!!')
catch
    WML = [];
    mkdir([R.path.rootn '\outputs\' R.out.tag ]);
    save([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\MultiStartListWML'],'WML')
    disp('Making Mod List!!')
end
N = 8; % number of starts
%% Prepare Model
for multiStart = 1:2*N
    load([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\MultiStartListWML'])
    if ~any(intersect(WML,multiStart))
        WML = [WML multiStart];
        save([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\MultiStartListWML'],'WML')
        disp('Writing to Mod List!!')
        fprintf('Now Fitting Mulitstart %.0f',multiStart)
        f = msgbox(sprintf('Fitting Multistart %.0f',multiStart));
        
        if multiStart < (N+1)
            R.data.feat_emp = feat_sim{1};
            R.data.feat_xscale{1} = R.frqz;
        else
            R.data.feat_emp = feat_sim{2};
            R.data.feat_xscale{1} = R.frqz;
        end
        %% Check Data
        R.plot.outFeatFx({R.data.feat_emp},{},R.data.feat_xscale,R,1,[])
        modelspec = eval(['@MS_rat_STN_GPe_ModComp_Model' num2str(1)]);
        [R,p,m] = modelspec(R);
        pause(2)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        R.out.dag = sprintf('NPD_STN_GPe_MultiStart_M%.0f',multiStart); % 'All Cross'
        % Delete Previous Saves
        delete([R.path.rootn 'outputs\' R.out.tag '\' R.out.dag '\*'])
        R.SimAn.rep = 512;
        R = setSimTime(R,32);
        R.plot.save = 0;
        R.Bcond = 0;
        R.plot.flag = 1;
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
