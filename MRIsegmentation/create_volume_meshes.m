function create_volume_meshes(input_img, maxvoxelvolume, input_coordsys, ...
                              output_coordsys)
%
% Input:        input_img <string> fullpath to folder, where the results of
%                                  mysegment.m are stored
%               maxvoxelvolume <float> desired maximal volume per tetrahedra
%               input_coordsys <string> name of coordinate system of input_img
%                                       (optional).
%               output_coordsys <string> name of desired coordinate system of
%                                        created mesh (optional).
%
% This function creates tetrahedral volume meshes from scalp, skull, csf, white
% matter and gray matter on the output of mysegment.m. It does the following:
% 1. Loads the tissue masks and combines gray and white matter to whitegray
% 2. Saves the final segmentation as segmentedmri
% 3. Calls iso2mesh for creating a volume mesh for each tissue class.
% 4. Stores the resulting volume mesh.
%
% Example:
% create_volume_meshes('/tmp/head/1/'); expects outputs in directory
% /tmp/head/1/ with input files named: mask_air.nii, mask_skin.nii,
% mask_bone.nii, mask_csf.nii, mask_gray.nii, mask_white.nii
%
% (c) Nils Harmening, May 2020
% Neurotechnology group, Technische Universität Berlin, Germany

[dirname, base_filename, ext] = fileparts(input_img);
if nargin > 1
  maxvoxelvol = maxvoxelvolume;
else
  maxvoxelvol = 1;
end

%% Create segmentation struct
mri = ft_read_mri(input_img);
if nargin > 2
  segmentedmri.coordsys = input_coordsys;
else
  segmentedmri.coordsys = 'acpc';
end
segmentedmri.transform = mri.transform;
segmentedmri.unit = mri.unit;
segmentedmri.dim = mri.dim;

%% Load andycorrected segmentations
c1 = ft_read_mri(fullfile(dirname,'mask_air.nii'));
c2 = ft_read_mri(fullfile(dirname,'mask_skin.nii'));
c3 = ft_read_mri(fullfile(dirname,'mask_bone.nii'));
c4 = ft_read_mri(fullfile(dirname,'mask_csf.nii'));
c5 = ft_read_mri(fullfile(dirname,'mask_gray.nii'));
c6 = ft_read_mri(fullfile(dirname,'mask_white.nii'));
segmentedmri.air = c1.anatomy;
segmentedmri.scalp = c2.anatomy;
segmentedmri.skull = c3.anatomy;
segmentedmri.csf = c4.anatomy;
segmentedmri.gray = c5.anatomy;
segmentedmri.white = c6.anatomy;
clear c1 c2 c3 c4 c5 c6



%% Build volume meshes
if ~numel(dir(fullfile(dirname,strcat('mesh6_maxvox',num2str(maxvoxelvol),...
                                      '.mat'))))
    tissue = {'air', 'scalp','skull','csf','gray','white'};
    mesh = prepare_mesh_iso2mesh(segmentedmri, tissue, maxvoxelvol);
    save(fullfile(dirname, strcat('mesh6_maxvox', num2str(maxvoxelvol),...
                                 '.mat')), 'mesh');
    savemsh(mesh.pos, [mesh.tet mesh.tetlabel], ...
            fullfile(dirname, strcat('mesh6_maxvox', num2str(maxvoxelvol), ...
                                     '.msh')), mesh.tissues);
end


%% Build volume meshes for HArtMuT
segmentedmri = prepend_terminal_slices(segmentedmri, 46, 5, 'zeros'); %'copyterminal');
save(fullfile(dirname,'segmentedmri_hartmut'), 'segmentedmri');

if ~numel(dir(fullfile(dirname,strcat('mesh6_maxvox',num2str(maxvoxelvol),...
                                      '_hartmut.mat'))))
    tissue = {'air', 'scalp','skull','csf','gray','white'};
    mesh = prepare_mesh_iso2mesh(segmentedmri, tissue, maxvoxelvol);
    save(fullfile(dirname, strcat('mesh6_maxvox', num2str(maxvoxelvol),...
                                 '_hartmut.mat')), 'mesh');
    savemsh(mesh.pos, [mesh.tet mesh.tetlabel], ...
            fullfile(dirname, strcat('mesh6_maxvox', num2str(maxvoxelvol), ...
                                     '_hartmut.msh')), mesh.tissues);
end
end %create_volume_meshes




function seg2 = prepend_terminal_slices(seg, nPad, nSlices, padModeOther)
% Like previous function, but in padModeOther='copyterminal' it ALSO fills the
% original inferior empty slices (1..kT-1) for each non-scalp tissue with its
% own terminal content (slice kT), so tissues don't "disappear" below.
%
% Assumes dim3 is inferior->superior increasing.

if nargin < 4 || isempty(padModeOther)
  padModeOther = 'zeros';
end

assert(isfield(seg,'scalp') && ndims(seg.scalp)==3, 'seg.scalp must be a 3D volume.');
assert(isfield(seg,'transform') && all(size(seg.transform)==[4 4]), 'seg.transform (4x4) is required.');

seg2 = seg;

dim = size(seg.scalp);
nx = dim(1); ny = dim(2); nz = dim(3);

% tissue-like fields
fn = fieldnames(seg);
tissues = {};
for i=1:numel(fn)
  f = fn{i};
  v = seg.(f);
  if (islogical(v) || isnumeric(v)) && ndims(v)==3 && all(size(v)==dim)
    tissues{end+1} = f; %#ok<AGROW>
  end
end

% -------- SCALP: walking cumulative OR over first nSlices non-empty region --------
scalpBin = (seg2.scalp ~= 0);
has = squeeze(any(any(scalpBin,1),2));
k0 = find(has, 1, 'first');
if isempty(k0), error('Scalp segmentation is empty.'); end
kEnd = min(nz, k0 + nSlices - 1);

scalpClass = class(seg2.scalp);
cum = (seg2.scalp(:,:,kEnd) ~= 0);   % start from kEnd

for k = (kEnd-1):-1:k0
  cum = cum | (seg2.scalp(:,:,k) ~= 0);
  seg2.scalp(:,:,k) = cast(cum, scalpClass);
end

% fill all slices below k0 (including any empty slices) with final cum
if k0 > 1
  seg2.scalp(:,:,1:k0-1) = repmat(cast(cum, scalpClass), 1, 1, k0-1);
end

scalpPadSlice = cast(cum, scalpClass);

% -------- OTHER TISSUES: optionally fill their inferior empty region too --------
if strcmpi(padModeOther, 'copyterminal')
  for it = 1:numel(tissues)
    t = tissues{it};
    if strcmpi(t,'scalp')
      continue;
    end
    vol = seg2.(t);

    bin = (vol ~= 0);
    hasT = squeeze(any(any(bin,1),2));
    kT = find(hasT, 1, 'first');

    if ~isempty(kT) && kT > 1
      % Fill 1..kT-1 with slice kT (this is the key correction)
      seg2.(t)(:,:,1:kT-1) = repmat(vol(:,:,kT), 1, 1, kT-1);
    end
  end
end

% -------- Prepend nPad slices to ALL tissues (keep alignment) --------
newNz = nz + nPad;

for it = 1:numel(tissues)
  t = tissues{it};
  vol = seg2.(t);
  cls = class(vol);

  if islogical(vol)
    volNew = false(nx, ny, newNz);
  else
    volNew = zeros(nx, ny, newNz, cls);
  end

  if strcmpi(t,'scalp')
    padSlice = scalpPadSlice;
  else
    switch lower(padModeOther)
      case 'zeros'
        if islogical(vol), padSlice = false(nx, ny);
        else,             padSlice = zeros(nx, ny, cls);
        end
      case 'copyterminal'
        % now that we've filled inferior region, slice 1 is guaranteed non-empty if tissue exists
        padSlice = vol(:,:,1);
      otherwise
        error('Unknown padModeOther: %s', padModeOther);
    end
  end

  volNew(:,:,1:nPad) = repmat(padSlice, 1, 1, nPad);
  volNew(:,:,nPad+1:end) = vol;

  seg2.(t) = volNew;
end

seg2.dim = [nx ny newNz];

% update transform once (we prepended for all tissues)
Shift = eye(4);
Shift(3,4) = -nPad;
seg2.transform = seg2.transform * Shift;

fprintf(['Scalp: walking OR fill k0=%d..kEnd=%d, filled below k0. ', ...
         'Other tissues: padModeOther=%s (and filled inferior empty region if copyterminal). ', ...
         'Prepended %d slices.\n'], ...
         k0, kEnd, padModeOther, nPad);

end

