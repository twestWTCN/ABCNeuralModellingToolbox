function modID = modelCompMaster_210320(Rorg,modlist,WML)
if nargin>2
    save([Rorg.rootn 'outputs\' Rorg.out.tag '\WorkingPermModList'],'WML')
end

%% Setup for parallelisation (multiple MATLAB sessions)
try
    load([Rorg.rootn 'outputs\' Rorg.out.tag '\WorkingPermModList'])
    disp('Loaded Perm Mod List!!')
    % If concatanating to previously computed model comp structure
catch
    WML = [];
    mkdir([Rorg.rootn 'outputs\' Rorg.out.tag ]);
    save([Rorg.rootn 'outputs\' Rorg.out.tag '\WorkingPermModList'],'WML')
    disp('Making Perm Mod List!!')
end
K = 0;
%% Main Loop
for modID = modlist
    K = K +1;
    load([Rorg.rootn 'outputs\' Rorg.out.tag '\WorkingPermModList'],'WML')
    permMod = [];
    if ~any(strmatch(modID,WML))
        WML = [WML modID];
        save([Rorg.rootn 'outputs\' Rorg.out.tag '\WorkingPermModList'],'WML')
        disp('Writing to PermMod List!!')
        fprintf(['Now Computing Probabilities for Model ' modID{1}])
        %         dagname = sprintf('NPD_InDrt_ModCompRed_M%.0f',modID);
        if Rorg.comptype == 1
            dagname = modID{1};
            SimMod = modID{1};
        elseif Rorg.comptype == 2
            SimData = Rorg.tmp.confmat(1,modID{1});
            SimMod = Rorg.tmp.confmat(2,modID{1});
            dagname = sprintf('NPD_STN_GPe_ConfMat_DataM%.0f_ParM%.0f',SimData,SimMod); % 'All Cross'
        end
        % Load Model
        load([Rorg.rootn 'outputs\' Rorg.out.tag '\' dagname '\modelspec_' Rorg.out.tag '_' dagname '.mat'])
        m = varo;
        % load modelfit
        load([Rorg.rootn 'outputs\' Rorg.out.tag '\' dagname '\modelfit_' Rorg.out.tag '_' dagname '.mat'])
        A = varo;
        p = A.BPfit;
        % Load Options
        load([Rorg.rootn 'outputs\' Rorg.out.tag '\' dagname '\R_' Rorg.out.tag '_' dagname '.mat'])
        R = varo;

        %         R.rootn = ['C:\Users\Tim\Documents\Work\GIT\SimAnneal_NeuroModel\Projects\' R.projectn '\'];
        %         R.rootm = 'C:\Users\Tim\Documents\Work\GIT\SimAnneal_NeuroModel\sim_machinery';
        
        R.Mfit = A;
%         R.Mfit.prior = prior;
        % load parbank?
        load([Rorg.rootn 'outputs\' Rorg.out.tag '\' dagname '\parBank_' Rorg.out.tag '_' dagname '.mat'])
        parBank =  varo;
        
        R.analysis.modEvi.eps = prctile(parBank(end,:),75);
        R.analysis.BAA.flag = 0; % Turn off BAA flag (time-locked analysis)
        parOptBank = parBank(1:end-1,parBank(end,:)>R.analysis.modEvi.eps);
        
        if  size(parOptBank,2)>1
            R.parOptBank = parOptBank;
            figure(K);
            R.analysis.modEvi.N = 1000;
            permMod = modelProbs(m.x,m,p,R);
        else
            permMod = [];
        end
        saveMkPath([Rorg.rootn 'outputs\' Rorg.out.tag '\' dagname '\modeProbs_' Rorg.out.tag '_' dagname '.mat'],permMod)
        pause(10)
        closeMessageBoxes
        close all
    end
end