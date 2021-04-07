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
R.path.datapath = [R.path.projpath '\data\Storage'];
% pathCell = regexp(path, pathsep, 'split'); onPath = any(strcmpi(spmpath, pathCell));
% if ~onPath; addpath(spmpath); spm eeg; close all; end

addpath(genpath([repopath '\ABC_dependencies']))
addpath(genpath([repopath '\sim_machinery']))
try; addpath(genpath([gitpath '\beta-burst-dyn'])); catch; warning('Optional Repo not found'); end
try; addpath(genpath([gitpath '\Spike-smr-reader'])); catch; warning('Optional Repo not found'); end
addpath(genpath(R.path.projpath))


