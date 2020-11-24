function pfldnm = getParFieldNames(p,m,prefix)
flag = 0;
if nargin<3
    prefix = '';
    flag = 1;
end
% Entry p is
X = p;
f   = fieldnames(X);
bCnt = 0;
pfldnm = {};
for fn = f'
    Y = X.(fn{1});
    fn{1}
    % If numeric
    if ~iscell(Y) && ~isstruct(Y)
        % If matrix
        if sum(size(Y)>1)>1
            for i = 1:size(Y,1)
                for j = 1:size(Y,2)
                    bCnt = bCnt + 1;
                    pfldnm{bCnt} = [fn{1} ' ' m.dipfit.model(i).source ' to ' m.dipfit.model(j).source];
                end
            end
        else % Simple vector
            for i = 1:numel(Y)
                bCnt = bCnt + 1;
                pfldnm{bCnt} = [prefix ' ' fn{1} num2str(i)];
            end
        end
    end
    
    % If structure then recurse
    if isstruct(Y)
        pfldn_rec = getParFieldNames(Y,m,fn{1});
        pfldnm = [pfldnm pfldn_rec];
        bCnt = numel(pfldnm);
    end
    
    % If there are subfields
    if iscell(Y)
        for L = 1:numel(Y)
            
            Z = Y{L};
            
            % Case of Matrix
            if sum(size(Z)>1)>1
                for i = 1:size(Z,1)
                    for j = 1:size(Z,2)
                        bCnt = bCnt + 1;
                        pfldnm{bCnt} = [fn{1} num2str(L) ' ' m.dipfit.model(i).source ' to ' m.dipfit.model(j).source];
                    end
                end
            end
            
            % Case of Structure
            if isstruct(Z)
                pfldn_rec = getParFieldNames(Z,m,[fn{1} ' ' m.dipfit.model(L).source]);
                pfldnm = [pfldnm pfldn_rec];
                bCnt = numel(pfldnm);
            end
            
        end
        
    end
end

if flag
    if numel(spm_vec(p))==numel(pfldnm)
        disp('Parameter name length matches!')
    else
        warning('Parameter name vector does not match parameter list!')
    end
end

