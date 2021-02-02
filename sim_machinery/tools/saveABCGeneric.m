function saveABCGeneric(R,varname,vartar)
for i = 1:numel(varname)
saveMkPath([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' R.out.dag '\' varname{i} '_' R.out.tag '_' R.out.dag '.mat'],vartar{i})
end