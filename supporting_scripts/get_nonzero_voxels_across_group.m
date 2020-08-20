function [masked_data,gmmask_composite] = get_nonzero_voxels_across_group(data,gmmask)
% for traveling subjects: mask w gm_mask, then calculate common mask of all nonzero voxels and re-mask
% gmmask = load_untouch_nii('/Users/stephanie/Documents/data/traveling_subs/5thpass_template_GM_WM_CSF_resl_crop_resampled.nii.gz');

data_gmmask=get_masked_data(data,gmmask);
data_gmmask_mat=cell2mat(data_gmmask);

for i=1:size(data_gmmask_mat,1)
row_zero_counts(i)=sum(data_gmmask_mat(i,:)==0);
% row_zero_counts(i)=abs(sum(data_gmmask_mat(i,:))<0.01);
end

datamask=row_zero_counts==0;
masked_data=get_masked_data(data_gmmask,datamask);

% return mask
gmmask_ids=find(gmmask);
gmmask_composite=gmmask;
for i=1:size(gmmask_ids)
    gmmask_composite(gmmask_ids(i))=gmmask(gmmask_ids(i))*datamask(i);
end

end

