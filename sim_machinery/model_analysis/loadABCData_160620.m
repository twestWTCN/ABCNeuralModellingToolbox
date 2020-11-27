function [Rout,m,p,parBank,permMod,parHist,bankSave,kldHist] = loadABCData_160620(R)
% Load Options
load([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\' R.out.dag '\R_' R.out.tag '_' R.out.dag '.mat'])
Rout = varo;
try; Rout.path.BBA_path = R.BBA_path; end
try; Rout.path.rootn = R.path.rootn; end
try;Rout.path.rootm = R.path.rootm; end
try;Rout.path.filepathn = R.path.filepathn; end

if nargout>1
    % Load Model
    load([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\' R.out.dag '\modelspec_' R.out.tag '_' R.out.dag '.mat'])
    m = varo;
end
if nargout>2
    % load modelfit
    load([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\' R.out.dag '\modelfit_' R.out.tag '_' R.out.dag '.mat'])
    Rout.Mfit = varo;
    p = Rout.Mfit.Pfit;
    
end
if nargout>3
    % Load Parbank
    load([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\' R.out.dag '\parBank_' R.out.tag '_' R.out.dag '.mat'])
    parBank =  varo;
end
if nargout>4
    try
        load([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' R.out.dag '\modeProbs_' R.out.tag '_'  R.out.dag '.mat'])
        permMod = varo; %i.e. permMod
    catch
        warning('ModelProbs was not available')
        permMod = [];
    end
end
if nargout>5
    load([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\' R.out.dag '\parHist_' R.out.tag '_' R.out.dag '.mat'])
    parHist = varo;
end
if nargout>6
    load([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\' R.out.dag '\bankSave_' R.out.tag '_' R.out.dag '.mat'])
    bankSave = varo;
end
if nargout>7
    load([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\' R.out.dag '\klHist_' R.out.tag '_' R.out.dag '.mat'])
    kldHist = varo;
end

