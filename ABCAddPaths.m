function R = ABCAddPaths(repopath,routname)
% This function sets up your path and points to the right toolboxes/etc
% restoredefaultpath
switch getenv('computername')
    case 'DESKTOP-94CEG1L'
        usname = 'timot';
        gitpath =  'C:\Users\timot\Documents\GitHub'; % points towards github repos
%         spmpath = 'C:\Users\timot\Documents\GitHub\spm12'; % points towrds
    case 'DESKTOP-4VATHIO'
        gitpath =  'D:\GITHUB'; % points towards github repos
%         spmpath = 'C:\Users\timot\Documents\GitHub\spm12'; % points towrds


    otherwise
        error('You need to setup your ABCAddPaths! Add paths for your computer...')
end
R.path.root = [repopath];
R.path.rootn = R.path.root; 
R.path.projectn = routname;
R.path.projpath =  [R.path.root '\Projects\' R.path.projectn];
R.path.datapath =  [R.path.projpath filesep 'data'];

addpath(genpath([gitpath '\ABCNeuralModellingToolbox\ABC_dependencies']))
addpath(genpath([gitpath '\ABCNeuralModellingToolbox\sim_machinery']))
addpath(genpath([repopath '\data']));
addpath(genpath([repopath '\model_fx']));
addpath(genpath([repopath '\ModelSpecs']));
addpath(genpath([repopath '\priors']));
addpath(genpath([repopath '\plotting']));
addpath(genpath([repopath '\routine\' routname]))
addpath(genpath([repopath '\statsfx']))
addpath(genpath([repopath '\dependencies']))
addpath(genpath(R.path.projpath))

try; addpath(genpath([gitpath '\beta-burst-dyn'])); catch; warning('Optional Repo not found'); end
try; addpath(genpath([gitpath '\Spike-smr-reader'])); catch; warning('Optional Repo not found'); end

