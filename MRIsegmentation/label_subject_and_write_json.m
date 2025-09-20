function label_subject_and_write_json(input_img, atlas_dir, num_vertices, out_json, opts)
% LABEL_SUBJECT_AND_WRITE_JSON
%
% Input:
%   input_img     <string>  Full path to a subject T1 (or any NIfTI in the CAT12 output
%                           folder). Its folder is used to locate CAT12 surfaces and its
%                           affine is used to transform the final OBJ to voxel space.
%   atlas_dir     <string>  Folder containing the precomputed fsaverage6 assets (created
%                           once via the Python helper):
%                               fsaverage6.L.sphere.reg.surf.gii
%                               fsaverage6.R.sphere.reg.surf.gii
%                               fsaverage6.L.Schaefer17_<N>.label.gii
%                               fsaverage6.R.Schaefer17_<N>.label.gii
%                               fsaverage6.L.Schaefer17_<N>.labelmap.json   (optional)
%                               fsaverage6.R.Schaefer17_<N>.labelmap.json   (optional)
%   num_vertices  <int>     (Optional) Target total vertex count for the combined cortex.
%                           If provided and smaller than the current total, each hemisphere
%                           is downsampled independently to floor(0.5 * num_vertices)
%                           **before** labeling and export. If omitted or larger than the
%                           current vertex count, no downsampling is performed.
%   out_json      <string>  (Optional) Output path for the parcel metadata JSON.
%                           Default: <folder_of_input_img>/parcels.json
%   opts          <struct>  (Optional) Behavior switches / overrides:
%       .ZeroIndexed        <bool>  Default true. If true, vertex indices in parcels.json
%                                   start at 0 (common for Python). If false, they start at 1.
%       .IncludeMedialWall  <bool>  Default true. Keep labels <= 0 (medial wall/background).
%       .ApplyInverseAffine <bool>  Default true. Apply inv(NIfTI affine) to surface vertices
%                                   before writing mask_brain.obj (OBJ ends up in voxel coords).
%       .NameMapL / .NameMapR       containers.Map(int32->string) to override label names per hemi.
%       .ColorMapL / .ColorMapR     containers.Map(int32->[r g b]) to override colors per hemi (0..1).
%
% This function takes CAT12’s subject-specific cortical surfaces and assigns Schaefer
% parcellation labels to every vertex via **spherical nearest-neighbor** transfer:
%   1) Reads CAT12 outputs from <folder_of_input_img>/surf:
%        - lh.central.<name>.gii, rh.central.<name>.gii
%        - lh.sphere.reg.<name>.gii, rh.sphere.reg.<name>.gii
%   2) (Optional) Downsamples each hemisphere (central + sphere) to match the requested
%      vertex budget while preserving 1:1 central↔sphere correspondence.
%   3) Transfers Schaefer labels from fsaverage6 to the subject by NN on the registered
%      sphere (fold-informed alignment), using the static fsaverage6 files in atlas_dir.
%   4) Concatenates LH then RH central surfaces, applies inv(NIfTI affine) from input_img,
%      and writes a single high-resolution cortex as:
%          <folder_of_input_img>/mask_brain.obj
%   5) Writes a pretty-printed JSON with one entry per parcel:
%          <folder_of_input_img>/parcels.json
%      Each entry contains:
%          - "Vertices": indices (into mask_brain.obj) of all parcel vertices
%          - "Color":    [r, g, b] in 0..1 (from labelmap or overrides; fallback gray)
%          - "Label":    parcel name (from labelmap or overrides; fallback "ID_<k> <hemi>")
%          - "Region":   "LH" or "RH"
%      Note: The optional *.labelmap.json files distributed with the atlas provide keys,
%      human-readable names (e.g., "Vis_1_LH"), and RGB colors. If missing, placeholders
%      are used and you may override with opts.*Map*.
%
% Requirements:
%   - SPM on path (for spm_vol to get the NIfTI affine).
%   - CAT12 surfaces present in <folder_of_input_img>/surf with the <name> derived from
%     input_img’s filename.
%   - fsaverage6 atlas assets in atlas_dir. Other parcellation atlasses can also be prepared.
%     Please contact me for this.
%
% Outputs (written next to input_img unless out_json is specified):
%   - mask_brain.obj        Combined cortex (LH then RH). OBJ uses 1-based faces (spec),
%                           vertices transformed to voxel coordinates if ApplyInverseAffine=true.
%   - parcels.json          Pretty-printed parcel metadata compatible with downstream tools.
%
% Example:
%   % No downsampling, zero-indexed JSON (default)
%   label_subject_and_write_json('/tmp/subj01/T1.nii', './Schaefer2018_Parcellations');
%   % Downsample to ~20k total vertices (≈10k per hemi) before labeling
%   label_subject_and_write_json('/tmp/subj01/T1.nii', './Schaefer2018_Parcellations', 20000);
%
% Notes:
%   - Vertex indices in parcels.json refer to mask_brain.obj with LH vertices first, then RH.
%   - The labeling is categorical NN on the sphere (no smoothing/mode filters applied).
%   - If you provide NameMap*/ColorMap* in opts, keys must match the atlas label keys
%     (typically FreeSurfer struct_id integers).
%
% This function is typically called after start_cat in the segmentation pipeline and
% produces the final surface + parcel metadata used for visualization and analysis in
% cedalions TwoSurfaceHeadModel.from_surfaces(...)
%
% (c) Nils Harmening, September 2025
% Neurotechnology group and IBS-lab, Technische Universität Berlin, Germany



[dirname, name, ~] = fileparts(input_img);
% --- optional: downsample each hemi BEFORE labeling ---
target_hemi = [];
if nargin >= 3 && ~isempty(num_vertices) && isnumeric(num_vertices) && isfinite(num_vertices)
    target_hemi = floor(0.5 * double(num_vertices));
end

if nargin < 4 || isempty(out_json), out_json  = fullfile(dirname,'parcels.json'); end
if nargin < 5, opts = struct; end
if ~isfield(opts,'ZeroIndexed'),       opts.ZeroIndexed = true; end
if ~isfield(opts,'IncludeMedialWall'), opts.IncludeMedialWall = true; end
if ~isfield(opts,'ApplyInverseAffine'),opts.ApplyInverseAffine = true; end

% --- load subject (CAT12) ---
gL = gifti(fullfile(dirname,'surf', ['lh.central.'    name '.gii']));
gR = gifti(fullfile(dirname,'surf', ['rh.central.'    name '.gii']));
sL = gifti(fullfile(dirname,'surf', ['lh.sphere.reg.' name '.gii']));
sR = gifti(fullfile(dirname,'surf', ['rh.sphere.reg.' name '.gii']));
assert(size(gL.vertices,1)==size(sL.vertices,1),'LH central/sphere mismatch');
assert(size(gR.vertices,1)==size(sR.vertices,1),'RH central/sphere mismatch');


if ~isempty(target_hemi) && target_hemi > 0
    % LH
    nL0 = size(gL.vertices,1);
    if target_hemi < nL0
        bndL.pos = double(gL.vertices); bndL.tri = double(gL.faces); bndL.Comment = 'LH';
        [NewL, I_L] = tess_downsize(bndL, target_hemi);
        gL.vertices = NewL.pos; gL.faces = NewL.tri;
        % keep the same vertices on the sphere for 1:1 mapping
        sL.vertices = sL.vertices(I_L, :);
        if isprop(sL,'faces') && ~isempty(sL.faces)
            % optional: reindex sphere faces if present (not needed for labeling)
            sL.faces = reindex_faces(sL.faces, I_L, size(NewL.pos,1));
        end
    end
    % RH
    nR0 = size(gR.vertices,1);
    if target_hemi < nR0
        bndR.pos = double(gR.vertices); bndR.tri = double(gR.faces); bndR.Comment = 'RH';
        [NewR, I_R] = tess_downsize(bndR, target_hemi);
        gR.vertices = NewR.pos; gR.faces = NewR.tri;
        sR.vertices = sR.vertices(I_R, :);
        if isprop(sR,'faces') && ~isempty(sR.faces)
            sR.faces = reindex_faces(sR.faces, I_R, size(NewR.pos,1));
        end
    end
end

% --- load fsaverage sphere + labels (static) ---
fsL = gifti(fullfile(atlas_dir,'fsaverage6.L.sphere.reg.surf.gii'));
fsR = gifti(fullfile(atlas_dir,'fsaverage6.R.sphere.reg.surf.gii'));
labL_g = gifti(fullfile(atlas_dir,'fsaverage6.L.Schaefer17_600.label.gii')); % adjust atlas if needed
labR_g = gifti(fullfile(atlas_dir,'fsaverage6.R.Schaefer17_600.label.gii'));
fsL_v  = fsL.vertices;   fsR_v  = fsR.vertices;
fsL_lab = int32(labL_g.cdata(:));
fsR_lab = int32(labR_g.cdata(:));
assert(size(fsL_v,1)==numel(fsL_lab),'LH fsavg sphere/label length mismatch');
assert(size(fsR_v,1)==numel(fsR_lab),'RH fsavg sphere/label length mismatch');

% --- NN on registered sphere: fsavg → subject (now at possibly downsampled density) ---
kdtL = KDTreeSearcher(fsL_v);   idxL = knnsearch(kdtL, sL.vertices);
kdtR = KDTreeSearcher(fsR_v);   idxR = knnsearch(kdtR, sR.vertices);
subjLabL = fsL_lab(idxL);
subjLabR = fsR_lab(idxR);

% --- combined cortex (LH+RH) ---
V_L = gL.vertices; F_L = gL.faces;
V_R = gR.vertices; F_R = gR.faces;
nL  = size(V_L,1);
V   = [V_L; V_R];
F   = [F_L; (F_R + nL)];

% --- apply inverse MRI affine to vertices (→ voxel space), then write OBJ ---
obj_path = fullfile(dirname, 'mask_brain.obj');
if opts.ApplyInverseAffine
    try
        V = apply_inverse_affine_from_nifti(input_img, V);
    catch ME
        warning('Inverse affine failed; writing OBJ in surface/world mm coords.\n%s', ME.message);
    end
end
write_obj(obj_path, V, F);

% --- load (optional) label maps for names/colors (per hemi) (STRICT) ---
[nameKeyL, colorKeyL, nameOrdL, colorOrdL] = load_label_maps_strict(atlas_dir, 'L');
[nameKeyR, colorKeyR, nameOrdR, colorOrdR] = load_label_maps_strict(atlas_dir, 'R');

% --- build parcels (indices refer to the COMBINED OBJ) ---
parcels = {};
push_hemi('L', subjLabL, 0,  nameKeyL, colorKeyL, nameOrdL, colorOrdL);
push_hemi('R', subjLabR, nL, nameKeyR, colorKeyR, nameOrdR, colorOrdR);

% --- pretty JSON ---
jsonTxt = '';
try
    jsonTxt = jsonencode(parcels, 'PrettyPrint', true);
catch
    jsonTxt = pretty_json(jsonencode(parcels));
end
fid = fopen(out_json,'w');  assert(fid>0, 'Cannot open JSON for writing');
fwrite(fid, jsonTxt, 'char');  fclose(fid);

fprintf('Wrote %d parcels to %s\n', numel(parcels), out_json);
fprintf('Wrote transformed cortex to %s (%d vertices, %d faces)\n', obj_path, size(V,1), size(F,1));

% ---------------- nested helpers ----------------
    function push_hemi(hemiChar, subjLab, offset, nameByKey, colorByKey, nameByOrd, colorByOrd)
        u = unique(subjLab(:));
        if ~opts.IncludeMedialWall, u(u<=0) = []; end
        u = sort(u(:));
        for ii = 1:numel(u)
            lbl = u(ii);
            if lbl==0 && ~opts.IncludeMedialWall, continue; end
            where = find(subjLab==lbl) + offset;   % global indices (LH then RH)
            if opts.ZeroIndexed, where = where - 1; end

            if hemiChar=='L'
                name = try_name(opts,'NameMapL', lbl, nameByKey, nameByOrd, sprintf('ID_%d L',lbl));
                col  = try_color_strict(opts,'ColorMapL', lbl, colorByKey, colorByOrd, [0.5 0.5 0.5]);
            else
                name = try_name(opts,'NameMapR', lbl, nameByKey, nameByOrd, sprintf('ID_%d R',lbl));
                col  = try_color_strict(opts,'ColorMapR', lbl, colorByKey, colorByOrd, [0.5 0.5 0.5]);
            end
            region = iff(hemiChar=='L','LH','RH');

            parcels{end+1,1} = struct( ...
                'Vertices', where(:)', ...
                'Color',    col(:)', ...
                'Label',    name, ...
                'Region',   region ...
            );
        end
    end
end % main

% ---------- affine / IO helpers ----------
function Vout = apply_inverse_affine_from_nifti(nifti_path, Vin)
    Vhdr = spm_vol(nifti_path);      % requires SPM on path
    A    = Vhdr.mat;                 % voxel -> mm
    VinH = [Vin, ones(size(Vin,1),1)];
    VoutH = (inv(A) * VinH')';
    Vout = VoutH(:,1:3);
end

function write_obj(path, V, F)
    fid = fopen(path,'w'); assert(fid>0, 'Cannot write OBJ');
    for i=1:size(V,1), fprintf(fid, 'v %.6f %.6f %.6f\n', V(i,1), V(i,2), V(i,3)); end
    for i=1:size(F,1), fprintf(fid, 'f %d %d %d\n', F(i,1), F(i,2), F(i,3)); end
    fclose(fid);
end

function Fnew = reindex_faces(Fold, keepIdx, nNew)
    % Map old vertex indices -> new [1..nNew] using keepIdx (sorted ascending)
    map = zeros(max(Fold(:)),1); map(keepIdx) = 1:numel(keepIdx);
    Fnew = map(Fold);
    if any(Fnew(:)==0), error('Face reindex failed: unseen vertex index.'); end
end

function s = pretty_json(compact)
    s = '';
    indent = 0; inStr = false; esc = false;
    ind = @() repmat(' ',1,indent*4);
    for ch = compact
        if inStr
            s(end+1) = ch; %#ok<AGROW>
            if esc, esc=false; elseif ch=='\', esc=true; elseif ch=='"', inStr=false; end
            continue
        end
        switch ch
            case '{', s=[s '{' newline ind()]; indent=indent+1;
            case '[', s=[s '[' newline ind()]; indent=indent+1;
            case '}', indent=indent-1; s=[s newline ind() '}'];
            case ']', indent=indent-1; s=[s newline ind() ']'];
            case ',', s=[s ',' newline ind()];
            case ':', s=[s ': '];
            case '"', s=[s '"']; inStr=true;
            otherwise, s=[s ch];
        end
    end
    s = string(s);
end

% ---------- label map loading (STRICT) ----------
function [nameByKey, colorByKey, nameByOrd, colorByOrd] = load_label_maps_strict(atlas_dir, hemiChar)
    nameByKey = []; colorByKey = []; nameByOrd = []; colorByOrd = [];
    mapfile = fullfile(atlas_dir, sprintf('fsaverage6.%s.Schaefer17_600.labelmap.json', hemiChar));
    if ~exist(mapfile,'file'), return; end
    M = jsondecode(fileread(mapfile));  % expects fields: key, name, color (0..1)
    n = numel(M);
    nameByOrd  = strings(n,1);
    colorByOrd = nan(n,3);
    for i = 1:n
        % name
        if ~isfield(M(i),'name') || isempty(M(i).name)
            error('Labelmap %s: entry %d missing "name".', mapfile, i);
        end
        nameByOrd(i) = string(M(i).name);
        % color
        if ~isfield(M(i),'color') || isempty(M(i).color)
            error('Labelmap %s: entry %d missing "color".', mapfile, i);
        end
        c = double(M(i).color);
        c = c(:).';
        if numel(c) ~= 3, error('Labelmap %s: entry %d color not 1x3.', mapfile, i); end
        if any(c < 0 | c > 1), error('Labelmap %s: entry %d color values must be in [0,1].', mapfile, i); end
        colorByOrd(i,:) = c;
    end
    % strong int32 key maps
    nameByKey  = containers.Map('KeyType','int32','ValueType','any');
    colorByKey = containers.Map('KeyType','int32','ValueType','any');
    for i=1:n
        if ~isfield(M(i),'key') || isempty(M(i).key)
            error('Labelmap %s: entry %d missing "key".', mapfile, i);
        end
        k = int32(M(i).key);
        nameByKey(k)  = nameByOrd(i);
        colorByKey(k) = colorByOrd(i,:);
    end
end

% ---------- tolerant map lookups + strict color ----------
function out = try_name(opts, mapField, key, nameByKey, nameByOrd, defaultVal)
    M = []; if isstruct(opts) && isfield(opts, mapField)
        tmp = opts.(mapField); if isa(tmp, 'containers.Map'), M = tmp; end
    end
    [hit, val] = map_lookup(M, key);
    if hit, out = string(val); return; end
    [hit, val] = map_lookup(nameByKey, key);
    if hit, out = string(val); return; end
    if ~isempty(nameByOrd) && double(key) >= 1 && double(key) <= numel(nameByOrd)
        out = nameByOrd(double(key)); return;
    end
    out = defaultVal;
end

function out = try_color_strict(opts, mapField, key, colorByKey, colorByOrd, defaultVal)
    M = []; if isstruct(opts) && isfield(opts, mapField)
        tmp = opts.(mapField); if isa(tmp, 'containers.Map'), M = tmp; end
    end
    [hit, val] = map_lookup(M, key);
    if ~hit, [hit, val] = map_lookup(colorByKey, key); end
    if ~hit && ~isempty(colorByOrd) && double(key) >= 1 && double(key) <= size(colorByOrd,1)
        val = colorByOrd(double(key),:); hit = true;
    end
    if ~hit, val = defaultVal; end
    val = double(val); val = val(:).';
    if numel(val) ~= 3, error('Color for label %d is not 1x3.', int32(key)); end
    if any(val < 0 | val > 1), error('Color for label %d out of [0,1].', int32(key)); end
    out = val;
end

function [hit, val] = map_lookup(M, key)
    hit = false; val = [];
    if isempty(M) || ~isa(M,'containers.Map'), return; end
    k_int = int32(key);
    if isKey(M, k_int), val = M(k_int); hit = true; return; end
    k_dbl = double(key);
    if isKey(M, k_dbl), val = M(k_dbl); hit = true; return; end
    K = keys(M);
    for i = 1:numel(K)
        Ki = K{i};
        if isnumeric(Ki) && double(Ki) == k_dbl
            val = M(Ki); hit = true; return;
        end
    end
end

% ---------- your requested downsampler (adapted from Brainstorm) ----------
function [NewTessMat, I, J] = tess_downsize(bnd, newNbVertices)
  % Downsize the new tesselation (Adapted from Brainstorm)
  NewTessMat = [];
  NewTessMat.Comment = bnd.Comment;
  I = [];
  J = [];
  oldNbVertices = size(bnd.pos, 1);
  if (newNbVertices >= oldNbVertices)
        disp(sprintf('TESS> Surface has %d vertices out of %d', oldNbVertices, newNbVertices));
        newNbVertices = oldNbVertices;
  end
  % Prepare variables
  bnd.tri = double(bnd.tri);
  bnd.pos = double(bnd.pos);
  dsFactor = newNbVertices / size(bnd.pos, 1);
  % Matlab's reducepatch: reduce number of vertices
  [NewTessMat.tri, NewTessMat.pos] = reducepatch(bnd.tri, bnd.pos, dsFactor);
  % Find the vertices that were kept by reducepatch
  [~, I, J] = intersect(bnd.pos, NewTessMat.pos, 'rows');
  % Re-order the vertices so that they are in the same order in the output surface
  [I, iSort] = sort(I);
  NewTessMat.pos = bnd.pos(I,:);
  J = J(iSort);
  % Re-order the vertices in the faces
  iSortFaces(J) = 1:length(J); %#ok<AGROW>
  NewTessMat.tri = iSortFaces(NewTessMat.tri);
  % Set the
  J = (1:length(J))';
end

function y = iff(cond, a, b), if cond, y = a; else, y = b; end, end
