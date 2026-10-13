function seg = restrict_to_head(seg, input_img, max_air_prior, tissues)
% RESTRICT_TO_HEAD  Drop tissue voxels that cannot be part of the head.
%
% SPM labels everything that is brighter than the (often extremely clean)
% background as tissue, including ghosting or averaging artefacts outside the
% head - e.g. the lateral smears of the ICBM152 template, which end up as
% floating blobs/fins in the scalp mesh. Two criteria remove them:
%
% 1. Population prior: voxels where the tissue probability map's air class
%    (sampled in native space through SPM's deformation field iy_<img>.nii)
%    exceeds max_air_prior are cleared. With the default 0.99 this only
%    touches locations where practically no head of the population reaches.
% 2. Connectivity: only the largest 26-connected component of the union of all
%    non-air tissues (the head) is kept.
%
% Input:   seg           <struct> binary tissue volumes as fields
%          input_img     <string> the segmented image; iy_<name>.nii and
%                                 <name>_seg8.mat (which records the TPM used)
%                                 next to it enable the prior criterion
%          max_air_prior <double> default 0.99; [] or >= 1 disables criterion 1
%          tissues       <cell>   fields to treat (default: those of scalp,
%                                 skull, csf, gray, white that exist)


if nargin < 3 || isempty(max_air_prior)
  max_air_prior = 0.99;
end
if nargin < 4
  tissues = {'scalp', 'skull', 'csf', 'gray', 'white'};
end
tissues = tissues(isfield(seg, tissues));
dim = size(seg.(tissues{1}));

%% 1. population prior
outside = false(dim);
if nargin >= 2 && ~isempty(input_img) && max_air_prior < 1
  [dirname, base, ext] = fileparts(input_img);
  % iy_<img>.nii is defined on the native grid and holds, per native voxel, the
  % corresponding MNI (template) coordinate in mm; y_<img>.nii is the reverse
  y_file = fullfile(dirname, ['iy_' base ext]);
  seg8 = fullfile(dirname, [base '_seg8.mat']);
  if exist(y_file, 'file') && exist(seg8, 'file')
    S = load(seg8);
    tpm_file = S.tpm(1).fname;
    if exist(tpm_file, 'file')
      Vair = spm_vol([tpm_file ',' num2str(numel(S.tpm))]);  % last class = air
      mni = zeros(3, prod(dim));
      for c = 1:3                                            % 5-D file: x, y, z components
        Vc = spm_vol(sprintf('%s,1,%d', y_file, c));
        mni(c, :) = reshape(spm_read_vols(Vc), 1, []);
      end
      vox = Vair.mat \ [mni; ones(1, size(mni, 2))];
      prior = spm_sample_vol(Vair, vox(1, :), vox(2, :), vox(3, :), 1);
      outside = reshape(prior > max_air_prior, dim);
      fprintf('restrict_to_head: air prior > %.2f (%s) at %d voxels\n', ...
              max_air_prior, tpm_file, nnz(outside));
    else
      warning('restrict_to_head:tpm', 'TPM %s not found, skipping the prior criterion', tpm_file);
    end
  else
    warning('restrict_to_head:deformation', ...
            'no %s / %s, skipping the prior criterion', y_file, seg8);
  end
end

%% 2. largest connected non-air component
union = false(dim);
for i = 1:numel(tissues)
  union = union | (seg.(tissues{i}) ~= 0 & ~outside);
end
cc = bwconncomp(union, 26);
n = cellfun(@numel, cc.PixelIdxList);
[~, k] = max(n);
head = false(dim);
head(cc.PixelIdxList{k}) = true;

for i = 1:numel(tissues)
  vol = seg.(tissues{i});
  removed = nnz(vol(~head));
  vol(~head) = 0;
  seg.(tissues{i}) = vol;
  if removed > 0
    fprintf('restrict_to_head: removed %d %s voxels outside the head\n', removed, tissues{i});
  end
end
fprintf('restrict_to_head: kept the largest of %d non-air components (%d voxels)\n', ...
        cc.NumObjects, n(k));
end
