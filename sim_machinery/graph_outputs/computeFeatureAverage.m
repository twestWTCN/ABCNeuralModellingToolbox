function feat_stat = computeFeatureAverage(R,feat_rep)
% this function compiles all of the model draws and finds statistics
% feat_stat is organized as [mean std sem iqr]
feat_stat = [];
for feat = 1:numel(R.data.datatype)
    Xtmp = [];
    % now collect terms of feature
    for rep = 1:numel(feat_rep)
        if ~isnan(feat_rep{rep}{1})
            if ndims(feat_rep{rep}{feat})>2
                Xtmp(:,:,:,:,:,rep) = feat_rep{rep}{feat}(1,1,1,R.siminds,:);
            else
                Xtmp(:,:,rep) = feat_rep{rep}{feat};
            end
            
        else % fill with nans
            if feat==1
                Xtmp(:,:,:,:,:,rep) = nan(1,1,1,1,size(R.data.feat_xscale{feat},2));
            else
                Xtmp(:,:,rep) = nan(1,size(R.data.feat_xscale{feat},2));
            end
        end
        
    end
    
    % now form average
    feat_stat{feat}(:,1) = squeeze(nanmean(Xtmp,ndims(Xtmp)));
    feat_stat{feat}(:,2) = squeeze(nanstd(Xtmp,[],ndims(Xtmp)));
    feat_stat{feat}(:,3) = squeeze(nanstd(Xtmp,[],ndims(Xtmp))./sqrt(numel(feat_rep)));
    feat_stat{feat}(:,4) = squeeze(iqr(Xtmp,ndims(Xtmp)));
    
end
