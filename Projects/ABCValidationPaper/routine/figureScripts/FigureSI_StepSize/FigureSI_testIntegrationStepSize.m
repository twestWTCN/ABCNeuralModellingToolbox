function FigureSI_testIntegrationStepSize(R)
%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE SII- Examination of simulatio step size
%%%%%%%%%%%%%%%%%%%%%%%%

clear ; close all
%This should link to your repo folder
repopath = 'C:\Users\Tim West\Documents\GitHub\ABCNeuralModellingToolbox';
addpath(repopath)
%This should be your projectname
projname = 'ABCValidationPaper';
R = ABCAddPaths(repopath,projname);

% Close all msgboxes
closeMessageBoxes

R.out.tag = 'figureSI_testIntegrationStepSize';
R.out.dag = '1';
fresh = 0;
if fresh
    r2 = zeros(1,4); feat_sim = cell(1,4);
    dts = [0.002 0.001 0.0005 0.0001];
    parfor dti = 1:4
        [r2(dti),pMAP,feat_sim{dti},xsims,xsims_gl{dti},wflag,Rmods(dti)] = getModelData(R,'figure5_ModelComp',10,1,dts(dti));
    end
    mkdir([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' R.out.dag])
    save([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' R.out.dag '\stepSizeInvestigation' R.out.tag '_' R.out.dag '.mat'],'xsims_gl','r2','feat_sim','Rmods')
else
    load([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' R.out.dag '\stepSizeInvestigation' R.out.tag '_' R.out.dag '.mat'],'xsims_gl','r2','feat_sim','Rmods')
end

cmap =  brewermap(6,'Reds');

for dti = 1:4
    Rmod = Rmods(dti);
    permMod.r2rep = r2(dti);
    permMod.feat_rep{1} = feat_sim{dti};
    h = figure(1);
    c = Rmod.data.feat_emp;
    Rmod.data.feat_emp = [];
    Rmod.data.feat_emp{1} = c;
    
    Rmod.plot.confint = 'none';
    PlotFeatureConfInt_gen170620(Rmod,permMod,h, cmap(dti+2,:));
    
end

figure(2)
histogram(r2)
xlabel('Pooled Mean Squared Error')
ylabel('Count')

figure(3)
cmap = linspecer(6);
for i = 1:25
    for j = 1:6
        plot(Rmod.IntP.tvec_obs,((xsims_gl{i}{1}(j,:))) - (12.*j),'Color',cmap(j,:))
        hold on
    end
end
xlabel('Time (s)')
ylabel('Channel')

figure(4)
cmap = linspecer(6);
for i = 1:25
    for j = 1:6
        plot(Rmod.IntP.tvec_obs,((xsims_gl{i}{1}(j,:))) - (12.*j),'Color',cmap(j,:))
        hold on
    end
end
xlabel('Time (s)')
ylabel('Channel')
xlim([30 35])