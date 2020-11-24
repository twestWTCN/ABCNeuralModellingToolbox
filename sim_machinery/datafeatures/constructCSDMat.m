function [F meancsd] = constructCSDMat(data,chloc_name,chlist,fsamp,N,R)
if isempty(N)
    N = floor(fsamp/R.obs.csd.df);
end
% N = R.obs.csd.pow2;
F_scale = R.frqz;
% Construct CSD matrix from data - take mean across channel replicates
for chloc = 1:size(chloc_name,2)
    chinds{chloc} = strmatch(chloc_name{chloc},chlist);
end
for chI = 1:size(chloc_name,2)
    for chJ = 1:size(chloc_name,2)
        for p = 1:size(chinds{chI},1)
            chindsP = chinds{chI};
            for r = 1:size(chinds{chJ},1)
                chindsR = chinds{chJ};
                [Pxy,F] = cpsd(data(chindsP(p),:),data(chindsR(r),:),hanning(2^N),[],2^N,fsamp);
                if nargin>5
                    Pxy = interp1(F,Pxy,F_scale);
                else
                    Pxy =  Pxy(F>4);
                end
                if istrue(R.obs.csd.ztranscsd)
                    Pxy = (Pxy-mean(Pxy))/std(Pxy);
                end
%                 if chI == chJ
                    if istrue(R.obs.csd.abovezero)
                        if min(Pxy)<0
                            Pxy = Pxy-min(Pxy);
                        else
                            Pxy = Pxy+min(Pxy);
                        end
                    end
%                 end
                xcsd(p,r,:) = Pxy;
            end
        end
        meancsd(chI,chJ,:) = squeeze(mean(mean(xcsd,1),2));
        clear xcsd
    end
end
F = F_scale;
% % take out symetrical CSD
% diaginds = [2 1; 3 1; 3 2];
% for i = 1:3
%     meancsd(diaginds(i,1),diaginds(i,2),:) = zeros(1,size(meancsd,3));
% end

end