function check_tissue_masks(dirname)
% CHECK_TISSUE_MASKS  Fail loudly when a segmentation has collapsed.
%
% SPM's unified segmentation can "succeed" while all brain classes end up
% empty (seen with the eTPM on the Colin27 template: everything skin/bone/air).
% mysegment then writes empty gray/white/csf masks and the failure only shows
% up much later as an obscure meshing error. Called by the Snakemake workflow
% right after mysegment.


tissues = {'gray', 'white', 'csf', 'bone', 'skin', 'air'};
n = zeros(1, numel(tissues));
for i = 1:numel(tissues)
  nii = load_untouch_nii(fullfile(dirname, ['mask_' tissues{i} '.nii']));
  n(i) = nnz(nii.img);
  total = numel(nii.img);
end
fprintf('tissue mask voxels: %s\n', ...
        strjoin(arrayfun(@(i) sprintf('%s=%d', tissues{i}, n(i)), 1:numel(tissues), ...
                         'UniformOutput', false), ', '));

if any(n(1:3) == 0) || n(5) > 0.5 * total
  error(['check_tissue_masks:collapsed', ...
         '\nThe tissue segmentation in %s is unusable (empty brain classes or ' ...
         'skin covering more than half of the volume). SPM most likely failed ' ...
         'on this image; check the input (intensities, orientation) or try ' ...
         'SPM''s standard TPM instead of the eTPM (config spm_tpm: TPM).'], dirname);
end
end
