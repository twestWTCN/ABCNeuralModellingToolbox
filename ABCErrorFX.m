function ABCErrorFX(R)
warning('Fitting failed')
try
    load([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\ERLIST'])
catch
    ERLIST{end+1} = [R.out.dag ' ' date ' ' lasterror.message];
    save([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\ERLIST'],'ERLIST')
end