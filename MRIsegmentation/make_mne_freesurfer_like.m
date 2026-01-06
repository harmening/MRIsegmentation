function make_mne_freesurfer_like(input_img, num_vertices)
%MAKE_MNE_FREESURFER_LIKE
%   Create an MNE/FreeSurfer-like folder structure from a segmented MRI:
%
%   Inputs:
%       input_img    : path to the original MRI (NIfTI)
%       num_vertices : integer used in the bnd4_<num_vertices>_corrected.mat
%
%   Assumes:
%     * All required files live in the same directory as input_img:
%         bnd4_<num_vertices>_corrected.mat    (FieldTrip new_bnd struct)
%       And CAT12 white-matter GIfTI surfaces:
%         surf/lh.white.<name>.gii
%         surf/rh.white.<name>.gii
%
%   Creates:
%       <dirname>/mne/<name>/mri/T1.nii
%       <dirname>/mne/<name>/bem/outer_skin.surf
%                                 outer_skull.surf
%                                 inner_skull.surf
%                                 inner_csf.surf  (extra; not required by MNE)
%       <dirname>/mne/<name>/surf/lh.white
%                                   rh.white
%
%   NOTE:
%       * No T1.mgz is created (pure MATLAB). Convert T1.nii to T1.mgz
%         yourself later with FreeSurfer or Python+nibabel if needed.
%       * No *-trans.fif is created here. You should generate it later
%         in MNE-Python via mne.gui.coregistration (or a script).

    %% Paths and folders
    [dirname, name, ~] = fileparts(input_img);

    subjects_dir = fullfile(dirname, 'mne');   % can be used as SUBJECTS_DIR
    subj_dir     = fullfile(subjects_dir, name);
    mri_dir      = fullfile(subj_dir, 'mri');
    bem_dir      = fullfile(subj_dir, 'bem');

    mkdir_if_needed(subjects_dir);
    mkdir_if_needed(subj_dir);
    mkdir_if_needed(mri_dir);
    mkdir_if_needed(bem_dir);

    %% 1) Copy NIfTI into mri/ as T1.nii (user can later convert to T1.mgz)
    if ~exist(input_img, 'file')
        error('Input NIfTI not found: %s', input_img);
    end
    T1_nifti_out = fullfile(mri_dir, 'T1.nii');
    fprintf('Copying MRI to %s\n', T1_nifti_out);
    copyfile(input_img, T1_nifti_out);

    %% 2) BEM surfaces from bnd4_<num_vertices>_corrected.mat (FieldTrip new_bnd)
    make_bem_surfaces_from_new_bnd(dirname, num_vertices, bem_dir);

    %% 3) Cortical white-matter surfaces from CAT12 GIfTI -> FreeSurfer .surf
    make_cortical_white_surfs(dirname, name, subj_dir);

    fprintf('Done.\n');
    fprintf('Subject folder: %s\n', subj_dir);
    fprintf(['NOTE: No T1.mgz and no *-trans.fif were created.\n' ...
             '      Please create them later in your Python/MNE pipeline.\n']);
end

%% ------------------------------------------------------------------------
function mkdir_if_needed(p)
    if ~exist(p, 'dir')
        mkdir(p);
    end
end

%% ------------------------------------------------------------------------
function make_bem_surfaces_from_new_bnd(dirname, num_vertices, bem_dir)
% Build BEM surfaces from FieldTrip new_bnd struct in mm and write
% FreeSurfer .surf files:
%   new_bnd(1) -> outer_skin.surf
%   new_bnd(2) -> outer_skull.surf
%   new_bnd(3) -> inner_skull.surf
%   new_bnd(4) -> inner_csf.surf (extra; cortex surface)

    matfile = fullfile(dirname, sprintf('bnd4_%d_corrected.mat', num_vertices));
    if ~exist(matfile, 'file')
        error('BEM .mat file not found: %s', matfile);
    end

    S = load(matfile, 'new_bnd');
    if ~isfield(S, 'new_bnd')
        error('File %s does not contain variable "new_bnd".', matfile);
    end
    new_bnd = S.new_bnd;

    % Convert units to mm using FieldTrip
    try
        new_bnd = ft_convert_units(new_bnd, 'mm');
    catch ME
        error(['ft_convert_units failed. Make sure FieldTrip is on the MATLAB path.\n' ...
               'Original error: %s'], ME.message);
    end

    % Mapping (following your description):
    %   1: scalp   -> outer_skin.surf
    %   2: skull   -> outer_skull.surf
    %   3: csf     -> inner_skull.surf
    %   4: cortex  -> inner_csf.surf (extra)
    surf_outer_skin  = fullfile(bem_dir, 'outer_skin.surf');
    surf_outer_skull = fullfile(bem_dir, 'outer_skull.surf');
    surf_inner_skull = fullfile(bem_dir, 'inner_skull.surf');
    surf_inner_csf   = fullfile(bem_dir, 'inner_csf.surf'); % extra / not required by MNE

    fprintf('Repairing and writing BEM FreeSurfer surfaces from new_bnd:\n');

    for ii = 1:4
        [new_bnd(ii).pos, new_bnd(ii).tri] = repair_bem_surface(new_bnd(ii).pos, new_bnd(ii).tri);
    end

    fprintf('Writing BEM FreeSurfer surfaces from new_bnd:\n');

    write_freesurfer_surf(surf_outer_skin,  new_bnd(1).pos, new_bnd(1).tri);
    write_freesurfer_surf(surf_outer_skull, new_bnd(2).pos, new_bnd(2).tri);
    write_freesurfer_surf(surf_inner_skull, new_bnd(3).pos, new_bnd(3).tri);
    write_freesurfer_surf(surf_inner_csf,   new_bnd(4).pos, new_bnd(4).tri);
end

%% ------------------------------------------------------------------------
function make_cortical_white_surfs(dirname, name, subj_dir)
% Convert CAT12 GIfTI white-matter and sphere surfaces to FreeSurfer .surf:
%   dirname/surf/lh.white.<name>.gii   -> subj_dir/surf/lh.white
%   dirname/surf/rh.white.<name>.gii   -> subj_dir/surf/rh.white
%   dirname/surf/lh.sphere.<name>.gii  -> subj_dir/surf/lh.sphere
%   dirname/surf/rh.sphere.<name>.gii  -> subj_dir/surf/rh.sphere

    surf_in_dir  = fullfile(dirname, 'surf');
    surf_out_dir = fullfile(subj_dir, 'surf');
    mkdir_if_needed(surf_out_dir);

    % ---- white surfaces ----
    lh_white_gii = fullfile(surf_in_dir, sprintf('lh.white.%s.gii',  name));
    rh_white_gii = fullfile(surf_in_dir, sprintf('rh.white.%s.gii',  name));

    if ~exist(lh_white_gii, 'file')
        warning('Left hemisphere white GIfTI not found: %s', lh_white_gii);
    else
        lh_white_out = fullfile(surf_out_dir, 'lh.white');
        gifti_to_freesurfer_surf(lh_white_gii, lh_white_out);
    end

    if ~exist(rh_white_gii, 'file')
        warning('Right hemisphere white GIfTI not found: %s', rh_white_gii);
    else
        rh_white_out = fullfile(surf_out_dir, 'rh.white');
        gifti_to_freesurfer_surf(rh_white_gii, rh_white_out);
    end

    % ---- sphere surfaces ----
    lh_sphere_gii = fullfile(surf_in_dir, sprintf('lh.sphere.%s.gii', name));
    rh_sphere_gii = fullfile(surf_in_dir, sprintf('rh.sphere.%s.gii', name));

    if ~exist(lh_sphere_gii, 'file')
        warning('Left hemisphere sphere GIfTI not found: %s', lh_sphere_gii);
    else
        lh_sphere_out = fullfile(surf_out_dir, 'lh.sphere');
        gifti_to_freesurfer_surf(lh_sphere_gii, lh_sphere_out);
    end

    if ~exist(rh_sphere_gii, 'file')
        warning('Right hemisphere sphere GIfTI not found: %s', rh_sphere_gii);
    else
        rh_sphere_out = fullfile(surf_out_dir, 'rh.sphere');
        gifti_to_freesurfer_surf(rh_sphere_gii, rh_sphere_out);
    end
end


%% ------------------------------------------------------------------------
function gifti_to_freesurfer_surf(gifti_file, surf_file)
% Read a GIfTI surface (CAT12 output) and write as FreeSurfer .surf.
% Expects a gifti object with .vertices and .faces, as in:
%   g = gifti(...);
%   g.vertices, g.faces

    if ~exist(gifti_file, 'file')
        error('GIfTI file not found: %s', gifti_file);
    end

    fprintf('Converting GIfTI -> FreeSurfer .surf:\n   %s\n   -> %s\n', ...
            gifti_file, surf_file);

    try
        g = gifti(gifti_file);
    catch ME
        error(['Failed to read GIfTI file with gifti(). ' ...
               'Make sure SPM / GIfTI toolbox is on the path.\n' ...
               'Original error: %s'], ME.message);
    end

    % --- directly use vertices & faces as you described --- %
    has_vertices = (isprop(g, 'vertices') || isfield(struct(g), 'vertices'));
    has_faces    = (isprop(g, 'faces')    || isfield(struct(g), 'faces'));

    if ~has_vertices || ~has_faces
        error(['Unexpected GIfTI structure for %s.\n' ...
               'I expected g.vertices and g.faces to exist.\n' ...
               'Manual test in MATLAB:\n' ...
               '   g = gifti(''%s'');\n' ...
               '   size(g.vertices), size(g.faces)\n'], gifti_file, gifti_file);
    end

    V = double(g.vertices);
    F = double(g.faces);

    write_freesurfer_surf(surf_file, V, F);
end

%% ------------------------------------------------------------------------
function write_freesurfer_surf(fname, verts, faces)
% Write a FreeSurfer .surf file (ASCII comments + binary data).
%   verts: N x 3 (mm, surface RAS / MRI coords)
%   faces: M x 3 (1-based indices)

    fprintf('Writing FreeSurfer surface: %s\n', fname);

    % Basic sanity checks
    if size(verts, 2) ~= 3
        error('verts must be N x 3, got %d x %d', size(verts,1), size(verts,2));
    end
    if size(faces, 2) ~= 3
        error('faces must be M x 3, got %d x %d', size(faces,1), size(faces,2));
    end

    % FreeSurfer uses 0-based vertex indices
    faces_zero_based = int32(faces - 1);

    fid = fopen(fname, 'wb', 'b'); % big-endian
    if fid == -1
        error('Cannot open output surface file for writing: %s', fname);
    end
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    % --- Magic number for *triangular* FreeSurfer surface: 0xFF 0xFF 0xFE --- %
    % This matches nibabel.freesurfer.io.write_geometry (magic = 16777214).
    fwrite(fid, [255 255 254], 'uint8');


    % Two comment lines (ASCII, newline-terminated)
    fprintf(fid, 'created by make_mne_freesurfer_like in MATLAB\n');
    fprintf(fid, '\n');

    nv = int32(size(verts, 1));
    nf = int32(size(faces_zero_based, 1));

    fwrite(fid, nv, 'int32');
    fwrite(fid, nf, 'int32');

    % Vertices: 3 x nv, float32 (stored as [x1...xn y1...yn z1...zn] OR [3 x nv])
    % nibabel/read_geometry expects FreeSurfer's standard: 3 x nv
    fwrite(fid, single(verts.'), 'float32');

    % Faces: 3 x nf, int32 (0-based)
    fwrite(fid, faces_zero_based.', 'int32');
end

%% ------------------------------------------------------------------------
function [V, F] = repair_bem_surface(V, F)
%REPAIR_BEM_SURFACE  Clean a BEM surface using iso2mesh.meshcheckrepair.
%
%   Fixes common issues:
%     - duplicated nodes/elements
%     - isolated nodes
%     - non-manifold vertices
%     - small holes / self-intersections (via meshfix if available)
%
%   Requires iso2mesh on the MATLAB path.

    if ~exist('meshcheckrepair', 'file')
        warning(['iso2mesh::meshcheckrepair not found on path, ' ...
                 'BEM surfaces are not repaired.']);
        return
    end

    % iso2mesh expects double
    V = double(V);
    F = double(F);

    % 1) Remove duplicates
    [V, F] = meshcheckrepair(V, F, 'dup');

    % 2) Remove isolated nodes
    [V, F] = meshcheckrepair(V, F, 'isolated');

    % 3) Remove non-manifold vertices (requires jmeshlib, part of iso2mesh)
    try
        [V, F] = meshcheckrepair(V, F, 'deep');
    catch
        warning('meshcheckrepair(...,''deep'') failed or jmeshlib missing; continuing.');
    end

    % 4) Try to repair/fill holes and self-intersections via meshfix
    try
        [V, F] = meshcheckrepair(V, F, 'meshfix');
    catch
        warning(['meshcheckrepair(...,''meshfix'') failed or meshfix binary missing.\n' ...
                 'Surface may still have holes; MNE can still complain.']);
    end

    % 5) Drop degenerate faces (same vertex index twice)
    bad = (F(:,1) == F(:,2)) | (F(:,1) == F(:,3)) | (F(:,2) == F(:,3));
    F(bad, :) = [];

    % 6) Optional: check for open boundaries (holes)
    if exist('surfedge', 'file')
        eg = surfedge(F);   % from iso2mesh
        if ~isempty(eg)
            warning('%d boundary edges remain in surface; mesh is still open.', size(eg, 1));
        end
    end
end

