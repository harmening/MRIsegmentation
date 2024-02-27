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

if nargin > 3
  segmentedmri.anatomy = mri.anatomy;
  segmentedmri = ft_convert_coordsys(segmentedmri, output_coordsys);
end

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
end %create_volume_meshes
