function [Rout,m,p,parBank] = loadABCData(R)
% Load Options
load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\R_' R.out.tag '_' R.out.dag '.mat'])
Rout = varo;
try; Rout.BBA_path = R.BBA_path; end
try; Rout.rootn = R.rootn; end
try;Rout.rootm = R.rootm; end
try;Rout.filepathn = R.filepathn; end

if nargout>1
    % Load Model
    load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelspec_' R.out.tag '_' R.out.dag '.mat'])
    m = varo;
end
if nargout>2
    % load modelfit
    load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelfit_' R.out.tag '_' R.out.dag '.mat'])
    Rout.Mfit = varo;
    p = Rout.Mfit.BPfit;
end
if nargout>3
    % Load Parbank
    load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\parBank_' R.out.tag '_' R.out.dag '.mat'])
    parBank =  varo;
end
% 'C:\Users\Tim\Documents\Work\GIT\SimAnneal_NeuroModel\Projects\Rat_NPD\outputs\InDrt_ModCompRev2\NPD_InDrt_ModCompRev2_M10\R_InDrt_ModCompRev2_NPD_InDrt_ModCompRev2_M10.mat'
% 'C:\Users\Tim\Documents\Work\GIT\SimAnneal_NeuroModel\Projects\Rat_NPD\outputs\InDrt_ModCompRev2\InDrt_ModCompRev2_M10\R_InDrt_ModCompRev2_InDrt_ModCompRev2_M10.mat'