clear ; close all;
%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 4- (II) Face Validation - Compute model evidence for each
%%%%%%%%%%%%%%%%%%%%%%%%

%This should link to your repo folder
repopath = 'C:\Users\timot\Documents\GitHub\ABCNeuralModellingToolbox';
cd(repopath)
addpath(repopath)
%This should be your projectname
projname = 'ABCValidationPaper';
R = ABCAddPaths(repopath,projname);

R.out.tag = 'figure4_2node_confusionMatrix';
R = ABCsetup_partII_FullModel(R);
closeMessageBoxes
fresh = 0;

%% Load Confusion Matrix Analysis
load([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\ConfData'],'feat_sim','pMAP','modlist'); %,'Rmod')
confmatlist = allcomb(modlist,modlist)';
R.tmp.confmat = confmatlist;

for i = 1:size(confmatlist,2)
    SimData = confmatlist(1,i);
    SimMod = confmatlist(2,i);
    daglist{i} = sprintf([R.out.tag '_DataM%.0f_ParM%.0f'],SimData,SimMod); % 'All Cross'
end
R.out.tag = 'figure2_FitDemo';
if fresh
    R.analysis.dagtype = 'arbitrary';
    modelCompMaster_160620(R,1:9,1:8,daglist) %,[1:8 10:12]
end
%% Analysis of model comp
PMOD = zeros(3,3);
DKL = zeros(3,3);
ACS = zeros(3,3);

iplist = {1:3,4:6,7:9};


for dataN = 1:3
    c = 1;
    ACCrep = {}; dkl = {}; SimData = []; SimMod = [];
    for modID = iplist{dataN}
        load([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' daglist{modID} '\modeProbs_' R.out.tag '_'  daglist{modID} '.mat'])
        permMod = varo; %i.e. permMod
        
            ACCrep{c} = permMod.ACCrep;
            dkl{c} = permMod.DKL;
        SimData(c) = confmatlist(1,modID);
        SimMod(c) = confmatlist(2,modID);
        c = c + 1;
    end
    
    % This gets the joint space epsilon that you can use to compute exceedence
    % probability
    prct = 50;
    ACCbankcat = horzcat(ACCrep{:});
    R.modcomp.modEvi.epspop = prctile(ACCbankcat,prct); % threshold becomes median of model fits
    
    
    pmod = cellfun(@(x) sum(x>=R.modcomp.modEvi.epspop)./(numel(x)+1),ACCrep,'UniformOutput',0);
    pmod = [pmod{:}];
    pModDist = (1-pmod)./sum(pmod);
    dkl = [dkl{:}];
    dklN = (numel(dkl)*dkl)./sum(dkl);
    
    PMOD(:,dataN) = pModDist'; %-log10(X)'; %X'; % the more positive the better
    CMP(:,dataN) = dklN'; %log10(dklN)'; %dkl'; the more negative the worse
    ACS(:,dataN) = (-log10(pModDist) - log10(dklN))'; % the more negative it is the better
    SM(:,dataN) = SimMod;
    SD(:,dataN) = SimData;
end

% cmap = brewermap(128,'*RdGy');

compStruc.pmod = PMOD;
compStruc.DKL = CMP;
compStruc.ACS = ACS;

cax = {[0 1],[0 2],[0 2]};
for n = 1:3
subplot(1,3,n)
h = plotModelComparisonMatrix(compStruc,n);
caxis(cax{n})
end
set(gcf,'Position',[625         613        1202         365])