function FigureSI_forwardUncertainty(R)
%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE SI- Examination of forward uncertainty of posterior model
%%%%%%%%%%%%%%%%%%%%%%%%
% Close all msgboxes
closeMessageBoxes

R.out.tag = 'figureSI_forwardUncertainty';
    R.out.dag = '1';
fresh = 0;
if fresh
    parfor i = 1:25
        [r2(i),pMAP,feat_sim{i},xsims,xsims_gl{i},wflag,Rmod{i}] = getModelData(R,'figure5_ModelComp',10);
    end
    mkdir([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' R.out.dag])
    save([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' R.out.dag '\multipleRealizaitions_' R.out.tag '_' R.out.dag '.mat'],'xsims_gl','r2','feat_sim','Rmod')
else
    load([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' R.out.dag '\multipleRealizaitions_' R.out.tag '_' R.out.dag '.mat'],'xsims_gl','r2','feat_sim','Rmod')
end

permMod.r2rep = r2;
permMod.feat_rep = feat_sim;
h = figure(1);
c = Rmod{1}.data.feat_emp;
Rmod{1}.data.feat_emp = [];
Rmod{1}.data.feat_emp{1} = c;

PlotFeatureConfInt_gen170620(Rmod{1},permMod,h, [1 0 0]);

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