function wflag = checkFileABCGeneric(R,varname)
for i = 1:numel(varname)
    if exist([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' R.out.dag '\' varname{i} '_' R.out.tag '_' R.out.dag '.mat'])
        wflag(i) = true;
    else
        wflag(i) = false;
    end
end