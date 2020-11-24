function pInd = parOptIndsPrec_110817(R,p,MN,set)
if nargin<4
    set =2;
end
plist = R.SimAn.pOptList;
for i = 1:length(plist)
    if ~isempty(strfind(plist{i},'src'))
        L = MN;
    else
        L = 1;
    end
    for src = 1:L
        if isequal(plist{i},'.A') || isequal(plist{i},'.B') % Case for matrices where we loop through
            for Ai = 1:size(eval(['p' plist{i}]),2)
                Y= eval(['p' plist{i} '{Ai}']);
                y = reshape(Y,1,[]);                
                
                X= eval(['p' plist{i} '_s{Ai}']);
                x = reshape(X,1,[]);
                xseq = rand(size(x));
                eval(['p' plist{i} '_s{Ai} = xseq;']);
                pvec = full(spm_vec(p));
                indK = [];
                for pin = 1:size(xseq,2)
                    indK(pin) = strfind(pvec',xseq(pin));
                end
                if set == 1
                    indK(y==-32) = NaN;
                elseif set == 2
                    indK(y==-32) = [];
                end
                eval(['pInd' plist{i} '_s{Ai} = indK;']);
                %                 eval(['pnew' plist{i} '{(Ai*2)-1} = X']);
            end
            
        else
            Y = eval(['p' plist{i}]);
            y = reshape(Y,1,[]);
            
            X = eval(['p' plist{i} '_s']);
            x = reshape(X,1,[]);
            xseq = rand(size(x));
            eval(['p' plist{i} '_s = xseq;']);
            pvec = full(spm_vec(p));
            indK = [];
            for pin = 1:size(xseq,2)
                indK(pin) = strfind(pvec',xseq(pin));
            end
            if set == 1
                indK(y==-32) = NaN;
            elseif set == 2
                indK(y==-32) = [];
            end
            eval(['pInd' plist{i} '_s = indK;']);
        end
    end
end




% vX = []
% f   = fieldnames(p);
% for i = 1:numel(f)
%     j = 1;
%     for j = 1:numel(p.(f{i}))
%         if iscell(p.(f{i}))
%             X = reshape(p.(f{i}){j},[],1);
%             for k = 0:size(X,1)-1
%             vX = [vX; X repmat(i,size(X)) repmat(j,size(X))];
%             vfield{i,j+k} = ['p.' f{i} '{' num2str(j) '}']
%             end
%         end
%         if isnumeric(p.(f{i}))
%             X = p.(f{i})(j)
%             vX = [vX; X i j];
%             vfield{i,j} = ['p.' f{i} '(' num2str(j) ')']
%         end
%         if isstruct(p.(f{i}))
%         end
%
%
%     end
%
% end
%     X = spm_vec({p.(f{i})});
%     vX = [vX; X repmat(i,size(X))];
% end
% pvec = full(spm_vec(p));
%
% for i = 1:length(plist)
%     ip = strmatch(plist{i}(2:3),f)
%         for j = 1:numel(p.(f{7}))
%
%
%
%