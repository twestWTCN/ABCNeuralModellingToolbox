function R = ABCForwardCompatibility(R)
flag = 0;
if ~isfield(R,'condN')
        R.condN = numel(R.condnames);
        flag = 1;
end
if isfield(R.IntP,'compFx')
    R.objfx.compFx = R.IntP.compFx;
    R.IntP = rmfield(R.IntP,'compFx');
    flag = 1;
end

if ~isfield(R.IntP,'fsample')
    R.IntP.fsample = 1/R.IntP.dt;
end

if flag
    warning('You appear to be using old/incomplete model specifications. Updated Fields!')
end