function varargout = loadABCGeneric(R,varname)
if ~iscell(varname)
    warning('Varnames must be given as cell array!')
    varname = {varname};
end
for i = 1:numel(varname)
    load([R.path.rootn '\outputs\' R.path.projectn '\'  R.out.tag '\' R.out.dag '\' varname{i} '_' R.out.tag '_' R.out.dag '.mat'])
    varargout{i} = varo;
end

