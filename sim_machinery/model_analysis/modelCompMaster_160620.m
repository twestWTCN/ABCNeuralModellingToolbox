function modID = modelCompMaster_160620(R,modlist,WML,daglist)
if nargin>2
    save([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\WorkingPermModList'],'WML')
end
closeMessageBoxes

if ~isfield(R.analysis,'dagtype')
    R.analysis.dagtype = 'normal';
end

% default to 1 realization of stochastic process
if ~isfield(R.SimAn,'RealzRep')
R.SimAn.RealzRep = 1;
end


%% Setup for parallelisation (multiple MATLAB sessions)
try
    load([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\WorkingPermModList'])
    disp('Loaded Perm Mod List!!')
    % If concatanating to previously computed model comp structure
catch
    WML = [];
    save([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\WorkingPermModList'],'WML')
    disp('Making Perm Mod List!!')
end

%% Main Loop
for modID = modlist
    load([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\WorkingPermModList'],'WML')
    permMod = [];
    if ~any(intersect(WML,modID))
        WML = [WML modID];
        save([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\WorkingPermModList'],'WML')
        disp('Writing to PermMod List!!')
        fprintf('Now Computing Probabilities for Model %.0f',modID)
        f = msgbox(sprintf('Probabilities for Model %.0f',modID));
        
        % Get Model Name
        switch R.analysis.dagtype
            case 'normal' % conventional model naming
                if isfield(R.out,'tag2')
                    dagcon = sprintf([R.out.tag2 '_M%.0f'],modID); % bugfix
                else
                    dagcon = sprintf([R.out.tag '_M%.0f'],modID);
                end
            case 'arbitrary' % used for confusion matrices
                dagcon = daglist{modID};
        end
        R.out.dag = dagcon;
        % Load Config
        varo = load([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' R.out.dag '\R_' R.out.tag '_' R.out.dag  '.mat'],'varo');
        
        % Replace with new version but maintain paths and tags
        tmp = varo.varo;
        tmp.path = R.path;
        tmp.plot = R.plot;
        tmp.out = R.out;
        tmp.analysis = R.analysis;
        %% Corrections to file structure to make compatible
        if ~iscell(tmp.data.feat_xscale)
            X = tmp.data.feat_xscale;
            tmp.data.feat_xscale = [];
            tmp.data.feat_xscale{1} = X;
        end
        if ~iscell(tmp.data.feat_emp)
            X = tmp.data.feat_emp;
            tmp.data.feat_emp = [];
            tmp.data.feat_emp{1} = X;
        end
        if ~iscell(tmp.data.datatype)
            X = tmp.data.datatype;
            tmp.data.datatype = [];
            tmp.data.datatype{1} = X;
        end
        
        if ~isfield(tmp,'chdat_name')
            tmp.chdat_name = tmp.chsim_name;
        end
        %%
        if ~isfield(R.objfx,'featweight')
            R.objfx.featweight = ones(size(R.data.datatype));
            warning('No feature weight specified so treating equally')
        end
        tmp.objfx.featweight = R.objfx.featweight;
        tmp.SimAn.RealzRep = R.SimAn.RealzRep;
        R  = tmp;
        
        [Rmod,m,p,parBank] = loadABCData_160620(R);
        if isfield(R.analysis,'comptype')
            if strcmp(R.analysis.comptype,'switch')
                Rmod.objfx.featweight = R.objfx.featweight;
                Rmod.SimAn.RealzRep = R.SimAn.RealzRep;
            end
        end
        % default to 1 realization of stochastic process
        if ~isfield(Rmod.SimAn,'RealzRep')
            Rmod.SimAn.RealzRep = 1;
        end
        R.analysis.modEvi.eps = parBank(end,R.SimAn.minRank);
        R.analysis.BAA.flag = 0; % Turn off BAA flag (time-locked analysis)
        R.out.dag = dagcon; %sprintf([R.out.tag '_M%.0f'],modID);
        permMod = modelProbs_160620(R,m.x,m,p,Rmod);
        saveMkPath([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' R.out.dag '\modeProbs_' R.out.tag '_' R.out.dag '.mat'],permMod)
        pause(10)
        closeMessageBoxes
        close all
    end
end