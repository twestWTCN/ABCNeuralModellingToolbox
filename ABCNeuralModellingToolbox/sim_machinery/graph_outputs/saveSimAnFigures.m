function saveSimAnFigures(R,ii)
saveallfiguresFIL_n([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' R.out.dag '\dataFeatures_' num2str(ii) '_'],'-jpg',1,'-r100',1);
saveallfiguresFIL_n([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' R.out.dag '\optimizationHistory_' num2str(ii) '_'],'-jpg',1,'-r100',2);
saveallfiguresFIL_n([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' R.out.dag '\posteriorEstimate_' num2str(ii) '_'],'-jpg',1,'-r100',3);