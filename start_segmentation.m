function start_segmentation(input_img, T2_optional)
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
%    triangular surface meshes for each tissue.
% 4. Based on the tissue masks, create_volume_meshes creates one tetrahedral
%    volume mesh including tissue labels for each mesh node, face, element.
% 5. Align_nonlinear nonlinearly aligns fiducials and electrodes of standard
%    mounting systems.
% 6. start_cat calls CAT12 for reliable gm and wm surface extraction.
% 7. CAT12's gm/wm output meshes are used by read_cat for building souremodels
%    of different size.
% 8. Electrodes, sourcemodels and meshes are transformed into ctf coordinate
%    system by transform_to_ctf.
%
% Example:
% segmentmri('/tmp/head/1/T1.nii'); 
% Segmentation results are stored at the location of the input_img and include 
% tissue masks ('mask...'), segmented MRI (segmentedmri.mat), 5-shell-surface-
% meshes (bnd5_1922.mat), extra cortex-smoothed 4-shell-surface-meshes
% (bnd4_1922.mat) and 5-tissue-volumetric mesh (mesh5_75000.mat), .....!!!!!...
% 
% (c) Nils Harmening, May 2020
% Neurotechnology group, Technische Universität Berlin, Germany

CWD = pwd; %erase(mfilename('fullpath'), 'start_segmentation'); %mfilename('fullpath');

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
num_vertices = 1922; % number of vertices of every surface mesh
%10000
maxvoxelvolume = 2; % max volume per tetrahedra of volume mesh
num_sources = 4000; % number of cortical sources in sourcemodel
%10000


%% Optional preprocessing + translation to ACPC (uncomment if you don't want to
%% preprocess)
if T2_optional
  preprocessing({input_img, T2_optional})
  [T2_filepath, T2_name, T2_ext] = fileparts(T2_optional);
  T2_optional = fullfile(T2_filepath, strcat(T2_name, '_RAS.nii'));
  input_img = fullfile(filepath, strcat(base_filename, '_RAS.nii'));
else
  T2_optional = [];
  preprocessing({input_img});
  input_img = fullfile(filepath, strcat(base_filename, '_RAS.nii'));
end



%% Start segmentation
Template = fullfile(CWD, 'Huang_et_al_2013', 'eTPM.nii');
normalize = false; 
start_seg(input_img, T2_optional, Template, normalize);

%% Andy's tools
[filepath, base_filename, ext] = fileparts(input_img);
isSmooth = true;
mysegment(filepath, base_filename, isSmooth);

%%% Create meshes
create_surface_meshes(input_img, num_vertices); 
%, input_coordsys, output_coordsys); 
create_volume_meshes(input_img, maxvoxelvolume);
%, input_coordsys, output_coordsys); 

%% Align electrodes and fiducials (nonlinear)
[dirname, base_filename, ext] = fileparts(input_img);
load(fullfile(dirname, strcat('bnd4_', num2str(num_vertices), ...
              '_corrected.mat')));
align_nonlinear(input_img, new_bnd(1));

%% Create sourcemodel
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
read_cat(cat12_path, input_img, num_sources);

%% Transform electrodes, sourcemodels and meshes to ctf coordinate system
transform_to_ctf(input_img, num_vertices, maxvoxelvolume, num_sources);

%% SPHARM
L = 50;
spharm(input_img, L);

end % start_segmentation
