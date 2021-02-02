function [flag,cnt] = spreadSession(R,cnt)
% This function can be used to split operation across multiple instances of
% MATLAB. Uses a shared working list WML.

flag = 0;

if nargin<2
    try
        load([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\WML'])
        disp('Loaded Mod List!!')
    catch
        WML = [];
        mkdir([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag ]);
        save([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\WML'],'WML')
        disp('Making Mod List!!')
    end
elseif isempty(cnt) || numel(cnt)>1 % this condition rewrites the WML
    WML = cnt;
    mkdir([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag ]);
    save([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\WML'],'WML')
    disp('Rewritting Mod List!!')    
elseif cnt>0
    load([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\WML'],'WML')
    if ~any(intersect(WML,cnt))
        WML = [WML cnt];
        save([R.path.rootn '\outputs\' R.path.projectn '\' R.out.tag '\WML'],'WML')
        disp('Writing to Mod List!!')
        f = msgbox(sprintf('Fitting Model to WML %.0f',cnt));
        flag = 1;
    end
elseif cnt == 0
    closeMessageBoxes;
    
end