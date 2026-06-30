function export_fieldtrip_hartmut_bnd(input_img, numvertices, outname)
% EXPORT_FIELDTRIP_HARTMUT_BND  Write the labelled bnd file the HArtMuT
% FieldTrip example (example_ft_dipolefitting_hartmut.m, cfg_source='mri') loads.
%
% Input:   input_img    <string> fullpath to the segmented T1 (same handle the
%                                 rest of the pipeline uses; only its folder is
%                                 read).
%          numvertices  <integer> surface vertex count used by the pipeline
%                                  (default 1922, must match start_segmentation).
%          outname      <string> output filename, default 'mri_bnd.mat',
%                                 written next to input_img.
%
% Writes one .mat holding standard FieldTrip constructs only:
%   bnd     4-element mesh array (scalp, skull, csf, cortex). Each element has
%           .pos (n-by-3 double), .tri (m-by-3 double), .unit 'mm',
%           .coordsys 'acpc'. The scalp is the neck-extended HArtMuT scalp.
%   tissue  {'scalp','skull','csf','cortex'}, aligned to bnd, the same cell array
%           ft_headmodel_openmeeg takes. The example keeps scalp/skull/csf and
%           uses csf as its inner boundary, so the cortex shell is shipped for
%           other uses only.
%   eye     optional. eye.pos is an n-by-3 set of ocular candidate positions for
%           ONE eye (left), pulled from the individualised HArtMuT artefact
%           sources, in the same acpc frame as bnd. Only written when the warped
%           artefact sourcemodel is present.
%
%
% Run after warp_HArtMuT_sources and before transform_to_ctf, so the surfaces
% and eye positions stay in acpc.
%
% Example:
%   export_fieldtrip_hartmut_bnd('/tmp/head/1/T1_RAS.nii', 1922);
%
% (c) Nils Harmening, June 2026
% Neurotechnology group, Technische Universität Berlin, Germany

if nargin < 2 || isempty(numvertices), numvertices = 1922;     end
if nargin < 3 || isempty(outname),     outname = 'mri_bnd.mat'; end

[dirname, ~, ~] = fileparts(input_img);
n = num2str(numvertices);

%% 1. Load the neck-extended surfaces
% Prefer the corrected (smoothed) hartmut surfaces, fall back to the raw ones.
candidates = { ...
  fullfile(dirname, ['bnd4_' n '_hartmut_corrected.mat']), ...
  fullfile(dirname, ['bnd4_' n '_hartmut.mat'])};
srcfile = '';
for i = 1:numel(candidates)
  if exist(candidates{i}, 'file') == 2, srcfile = candidates{i}; break; end
end
assert(~isempty(srcfile), ['No neck-extended hartmut surfaces found for ' n ...
  ' vertices. Run create_surface_meshes first.']);
S = load(srcfile);

% The corrected file stores 'new_bnd', the raw one stores 'bnd'.
if isfield(S, 'new_bnd')
  raw = S.new_bnd;
elseif isfield(S, 'bnd')
  raw = S.bnd;
else
  fn = fieldnames(S);
  raw = S.(fn{1});
end
assert(numel(raw) >= 4, 'Expected a 4-shell mesh array (scalp, skull, csf, cortex).');

%% 2. Normalise every shell to .pos / .tri / .unit / .coordsys
bnd = struct('pos', {}, 'tri', {}, 'unit', {}, 'coordsys', {});
for i = 1:4
  el = raw(i);
  if isfield(el, 'pos') && ~isempty(el.pos)
    p = el.pos;
  elseif isfield(el, 'pnt') && ~isempty(el.pnt)   % older FieldTrip meshes
    p = el.pnt;
  else
    error('Shell %d has neither .pos nor .pnt.', i);
  end
  bnd(i).pos      = double(p);
  bnd(i).tri      = double(el.tri);
  bnd(i).unit     = 'mm';
  bnd(i).coordsys = 'acpc';
end

%% 3. Tissue labels aligned to bnd
% Naming convention: 'brain' = 'csf' U 'cortex', so the innermost shell ships
% as 'cortex', not 'brain'.
tissue = {'scalp', 'skull', 'csf', 'cortex'};

%% 4. One eye's ocular candidates from the individualised artefact sources
eye = pull_eye_positions(dirname);  %#ok<NASGU>  may be empty

%% 5. Save
outpath = fullfile(dirname, outname);
if isempty(eye.pos)
  warning(['No warped artefact sourcemodel found, writing bnd and tissue ' ...
           'without eye.pos (the example treats eye as optional).']);
  save(outpath, 'bnd', 'tissue');
else
  save(outpath, 'bnd', 'tissue', 'eye');
end
fprintf('Wrote %s (scalp %d verts, eye candidates %d).\n', ...
  outpath, size(bnd(1).pos,1), size(eye.pos,1));

end % export_fieldtrip_hartmut_bnd


% ======================================================================
function eye = pull_eye_positions(dirname)
% Pick the left eye's ocular candidate positions from the warped HArtMuT
% artefact sourcemodel. The example mirrors this single eye onto the other side.
eye = struct('pos', zeros(0,3));

srcfile = fullfile(dirname, 'artefact_sourcemodel_HArtMuT_small.mat');
if exist(srcfile, 'file') ~= 2
  return;
end
A = load(srcfile);
if isfield(A, 'artefactmodel')
  am = A.artefactmodel;
elseif isfield(A, 'individual_artefactmodel')
  am = A.individual_artefactmodel;
else
  return;
end
if ~isfield(am, 'pos') || ~isfield(am, 'labels'), return; end

labels = am.labels;
if ~iscell(labels), labels = cellstr(labels); end

% Eye sources are labelled Eye... and carry a side: 'left', 'right' or the
% midline-symmetric 'leftright'. Keep one eye: 'left' but not 'leftright'.
isEye  = ~cellfun(@isempty, regexpi(labels, '^Eye', 'once'));
isLeft = ~cellfun(@isempty, regexpi(labels, 'left', 'once'));
isLR   = ~cellfun(@isempty, regexpi(labels, 'leftright', 'once'));
pick   = isEye & isLeft & ~isLR;

eye.pos = double(am.pos(pick, :));
end
