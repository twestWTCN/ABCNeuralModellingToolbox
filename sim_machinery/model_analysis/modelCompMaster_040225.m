function modID = modelCompMaster_040225(R,modlist,WML,daglist)
Rorg = R;
if nargin>2
    save([R.path.rootn '\outputs\' Rorg.path.projectn '\'  Rorg.out.tag '\WorkingPermModList'],'WML')
end
closeMessageBoxes

if ~isfield(Rorg.analysis,'dagtype')
    Rorg.analysis.dagtype = 'normal';
end

% default to 1 realization of stochastic process
if ~isfield(Rorg.SimAn,'RealzRep')
Rorg.SimAn.RealzRep = 1;
end

%% Setup for parallelisation (multiple MATLAB sessions)
try
    load([Rorg.path.rootn '\outputs\' Rorg.path.projectn '\'  Rorg.out.tag '\WorkingPermModList'])
    disp('Loaded Perm Mod List!!')
    % If concatanating to previously computed model comp structure
catch
    WML = [];
    save([Rorg.path.rootn '\outputs\' Rorg.path.projectn '\'  Rorg.out.tag '\WorkingPermModList'],'WML')
    disp('Making Perm Mod List!!')
end

%% Main Loop
for modID = modlist
    load([Rorg.path.rootn '\outputs\' Rorg.path.projectn '\'  Rorg.out.tag '\WorkingPermModList'],'WML')
    permMod = [];
    if ~any(intersect(WML,modID))
        WML = [WML modID];
        save([Rorg.path.rootn '\outputs\' Rorg.path.projectn '\'  Rorg.out.tag '\WorkingPermModList'],'WML')
        disp('Writing to PermMod List!!')
        fprintf('Now Computing Probabilities for Model %.0f',modID)
        f = msgbox(sprintf('Probabilities for Model %.0f',modID));
        
        % Get Model Name
        switch R.analysis.dagtype
            case 'normal' % conventional model naming
                if ~isfield(R.out,'tagexpression')
                    if isfield(R.out,'tag2')
                        dagcon = sprintf([Rorg.out.tag2 '_M%.0f'],modID); % bugfix
                    else
                        dagcon = sprintf([Rorg.out.tag '_M%.0f'],modID);
                    end
                else
                    dagcon = sprintf([Rorg.out.tagexpression],modID);
                end
            case 'arbitrary' % used for confusion matrices
                dagcon = daglist{modID};
        end
        R.out.dag = dagcon;
        %%
        if ~isfield(R.objfx,'featweight')
            R.objfx.featweight = ones(size(R.data.datatype));
            warning('No feature weight specified so treating equally')
        end
        
        [Rmod,m,p,parBank] = loadABCData_160620(R);
        if isfield(R.analysis,'comptype')
            if strcmp(R.analysis.comptype,'switch')
                if isfield(Rorg.objfx,'featweight')
                    Rmod.objfx.featweight = Rorg.objfx.featweight;
                end
                if isfield(Rorg.objfx,'errorFx')
                    Rmod.objfx.errorFx = Rorg.objfx.errorFx;
                end
                if isfield(Rorg.objfx,'RealzRep')
                    Rmod.objfx.RealzRep = Rorg.objfx.RealzRep;
                end
            end
        end
        % default to 1 realization of stochastic process
        if ~isfield(Rmod.SimAn,'RealzRep')
            Rmod.SimAn.RealzRep = 1;
        end
        Rmod.analysis.modEvi.N = R.analysis.modEvi.N;
        Rmod.analysis.modEvi.eps = parBank(end,Rmod.SimAn.minRank);
        Rmod.analysis.BAA.flag = 0; % Turn off BAA flag (time-locked analysis)
        permMod = modelProbs_160620(Rmod,m.x,m,p,Rmod);
        saveMkPath([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' R.out.dag '\modeProbs_' R.out.tag '_' R.out.dag '.mat'],permMod)
        % pause(10)
        closeMessageBoxes
        close all
    end
end