function start_segmentation_cedalion(input_img, T2_optional)
%
% Input:        input_img <string> fullpath to T1 MRI image to be segmented
%               T2_optional <string> fullpath to T2 weighted MRI image.
%                                    Supports SPM during segmentation.
%
% This function calls step-by-step the full segmentation pipeline (without)
% preprocessing. It does the following:
% 1. Start new_segment function of SPM12 with the eTPM.nii.
% 2. Andy's tools do an automatic clean up of the SPM output and store a mask
%    for each tissue.
% 3. Based on the tissue masks, create_surface_meshes creates and repaires 
%    triangular surface meshes for each tissue using iso2mesh and CGAL.
% 4. Based on the tissue masks, create_volume_meshes creates one tetrahedral
%    volume mesh including tissue labels for each mesh node, face, element
%    using iso2mesh and CGAL.
% 5. start_cat calls CAT12 for reliable gm and wm surface extraction.
% 6. CAT12's gm/wm output meshes are used by label_subject_and_write_json for
%    exporting a high-resolution cortex surface and adding SchaeferParcels to
%    it.
%
% Example:
% segmentmri('/tmp/head/1/T1.nii'); 
% Segmentation results are stored at the location of the input_img and include 
% tissue masks ('mask...'), segmented MRI (segmentedmri.mat), SPM scalp mesh
% are store as mask_scalp<num_verts>.stl and mask_scalp<num_verts>_smoothed.stl
% and SPM cortex mesh as mask_cortex<num_verts>.stl.
% CAT12 cortex surface is stored as mask_brain.obj and it's parcellation
% information as parcels.json.
% Nirfast meshes (mesh6_maxvox2.msh) in RAS for photon simulation using NIRFAST.
% 
% (c) Nils Harmening, September 2025
% Neurotechnology group and IBS-lab, Technische Universität Berlin, Germany

CWD = pwd;

% Specify Input
if nargin <2 || isempty(T2_optional) %no T2 specified
    T2_optional = [];
end
%input_img = fullfile(CWD, 'data', 'example.img');
disp(T2_optional);
%T2_optional = [];
[filepath, base_filename, ext] = fileparts(input_img);
if strcmp(ext, '.gz')
    gunzip(input_img, filepath);
    input_img = fullfile(filepath, base_filename);
    [filepath, base_filename, ext] = fileparts(input_img);
end

% DEFAULT MESH SIZE SETTINGS:
num_vertices = 15000; % number of vertices for scalp and SPM cortex mesh
num_cortex_verts = 15000; % number of vertices for CAT12 cortex mesh (and parcels)
maxvoxelvolume = 2; % max volume per tetrahedra of volume mesh (for NIRFASTER)


%% Start segmentation
Template = fullfile(CWD, 'Huang_et_al_2013', 'eTPM.nii');
normalize = false; 
start_seg(input_img, T2_optional, Template, normalize);

%% Andy's tools (segmentation postprocessing)
[filepath, base_filename, ext] = fileparts(input_img);
isSmooth = true;
mysegment(filepath, base_filename, isSmooth);

%% Create meshes
% for 2SHM
create_surface_meshes_cedalion(input_img, num_vertices); 
% for NIRFASTER
create_volume_meshes(input_img, maxvoxelvolume);


%% Create CAT12 cortex surface mesh with parcellation labeling (optional)
cat12_path = fullfile(CWD, 'fieldtrip', 'external', 'spm12', 'toolbox', ...
                      'cat12');
atlas = fullfile(cat12_path, 'atlases_surfaces', ...
                 'lh.aparc_DK40.freesurfer.annot');
%dartelTpm = fullfile(cat12_path, 'templates_1.50mm', ...
%                     'Template_1_IXI555_MNI152.nii'); % version 12
%dartelTpm = fullfile(cat12_path, 'templates_volumes',  ...
%                     'Template_1_IXI555_MNI152.nii'); % version 12.7
dartelTpm = fullfile(cat12_path, 'templates_MNI152NLin2009cAsym',  ...
                     'Template_1_Dartel.nii'); % version 12.8
start_cat(input_img, atlas, dartelTpm);

atlas_dir = fullfile(CWD, 'Schaefer2018_Parcellations');
label_subject_and_write_json(input_img, atlas_dir, num_cortex_verts);
end % start_segmentation
