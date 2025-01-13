function transform_to_ctf(input_img, numverts, maxvoxelvol, num_sources)
%
% Input:        input_img <string> fullpath to T1 MRI image to be segmented
%               numvertices <integer> number of desired mesh vertices.
%               maxvoxelvol <float> desired maximal volume per tetrahedra
%               num_sources <integer> number of sourcemodel points.
% 
% This function transforms all coordinate system depending results of this
% segmentaion pipeline into the ctf coordinate system. It does the follwoing:
% 1. Reads fiducials from nonlinear alignement.
% 2. Calculates the transformation.
% 3. Applies the transformation to surface and volume meshes, electrodes and
%    sourcemodels.
%
% Example:
% transform_to_ctf('/tmp/head/1/T1.nii');
% 
% (c) Nils Harmening, May 2020
% Neurotechnology group, Technische Universität Berlin, Germany

[dirname, base_filename, ext] = fileparts(input_img);
mkdir(fullfile(dirname, 'ctf'));

%% Load fiducials
load(fullfile(dirname, 'fiducials_nonlinear.mat'));
nas = fiducials_aligned.chanpos(1,:);
lpa = fiducials_aligned.chanpos(2,:);
rpa = fiducials_aligned.chanpos(3,:);

transform = ft_headcoordinates(nas, lpa, rpa, 'ctf');
save(fullfile(dirname, 'ctf', 'acpc2ctf.mat'), 'transform');

%% Convert meshes
cfg = [];
cfg.method        = 'fiducial';
cfg.coordsys      = 'ctf';
cfg.fiducial.nas  = nas; %position of NAS
cfg.fiducial.lpa  = lpa; %position of LPA
cfg.fiducial.rpa  = rpa; %position of RPA

new_bnd = [];
new_bnd.coordsys = 'ctf';
if numel(dir(fullfile(dirname, strcat('bnd4_', num2str(numverts), ...
                                      '_corrected.mat'))))
  load(fullfile(dirname, strcat('bnd4_', num2str(numverts), '_corrected.mat')));
  shells = {'scalp', 'skull', 'csf', 'cortex'};
  for i=1:4
    new_bnd_i = ft_meshrealign(cfg, new_bnd(i));
    new_bnd(i).pos = new_bnd_i.pos;
    new_bnd(i).tri = new_bnd_i.tri;
    om_save_tri(fullfile(dirname, 'ctf', strcat('bnd4_', num2str(numverts), ...
                                                '_corrected_', ...
                                                char(shells(i)), '.tri')), ...
                new_bnd(i).pos, new_bnd(i).tri);
  end
end
save(fullfile(dirname, 'ctf', strcat('bnd4_', num2str(numverts), ...
                                     '_corrected.mat')), 'new_bnd');

if numel(dir(fullfile(dirname, strcat('mesh6_maxvox', num2str(maxvoxelvol), ...
                                      '.mat'))))
  load(fullfile(dirname, strcat('mesh6_maxvox', num2str(maxvoxelvol), '.mat')));
  transform = ft_headcoordinates(nas, lpa, rpa, 'ctf');
  mesh = ft_transform_geometry(transform, mesh);
  mesh = rmfield(mesh, 'coordsys');
  % same as 
  %mesh = ft_meshrealign(cfg, mesh);
  save(fullfile(dirname, 'ctf', strcat('mesh6_maxvox', ...
                                       num2str(maxvoxelvol), '.mat')), 'mesh');
end
if numel(dir(fullfile(dirname, strcat('mesh5_maxvox', num2str(maxvoxelvol), ...
                                      '.mat'))))
  load(fullfile(dirname, strcat('mesh5_maxvox', num2str(maxvoxelvol), '.mat')));
  transform = ft_headcoordinates(nas, lpa, rpa, 'ctf');
  mesh = ft_transform_geometry(transform, mesh);
  mesh = rmfield(mesh, 'coordsys');
  % same as 
  %mesh = ft_meshrealign(cfg, mesh);
  save(fullfile(dirname, 'ctf', strcat('mesh5_maxvox', ...
                                       num2str(maxvoxelvol), '.mat')), 'mesh');
end

%% Convert cortex surface for SPHARM
transform = ft_headcoordinates(nas, lpa, rpa, 'ctf');

load(fullfile(dirname, 'surf', 'lh_surf.mat'));
bnd_lh = ft_transform_geometry(transform, bnd_lh);
save(fullfile(dirname, 'ctf', 'lh_surf.mat'), 'bnd_lh');
om_save_tri(fullfile(dirname, 'ctf', 'lh_surf.tri'), bnd_lh.pos, bnd_lh.tri);
load(fullfile(dirname, 'surf', 'lh_sphere.mat'));
lh_sphere = ft_transform_geometry(transform, lh_sphere);
save(fullfile(dirname, 'ctf', 'lh_sphere.mat'), 'lh_sphere');

load(fullfile(dirname, 'surf', 'rh_surf.mat'));
bnd_rh = ft_transform_geometry(transform, bnd_rh);
save(fullfile(dirname, 'ctf', 'rh_surf.mat'), 'bnd_rh');
om_save_tri(fullfile(dirname, 'ctf', 'rh_surf.tri'), bnd_rh.pos, bnd_rh.tri);
load(fullfile(dirname, 'surf', 'rh_sphere.mat'));
rh_sphere = ft_transform_geometry(transform, rh_sphere);
save(fullfile(dirname, 'ctf', 'rh_sphere.mat'), 'rh_sphere');



%% Convert electrodes
transform = ft_headcoordinates(nas, lpa, rpa, 'ctf');
elecfiles = {'standard_1005' 'standard_1010' 'standard_1020'};
for elecfile = elecfiles
  if numel(dir(fullfile(dirname, strcat('elec_', char(elecfile), '.mat'))))
    load(fullfile(dirname, strcat('elec_', char(elecfile), '.mat')));
    elec_aligned = ft_transform_geometry(transform, elec_aligned);
    save(fullfile(dirname, 'ctf', strcat('elec_', char(elecfile), '.mat')), ...
         'elec_aligned');
  end
end

%% Convert sourcemodels
transform = ft_headcoordinates(nas, lpa, rpa, 'ctf');
sourcefiles = {strcat('sourcemodel', num2str(num_sources))};
for sourcefile = sourcefiles
	load(fullfile(dirname, strcat(char(sourcefile), '.mat')));
	sourcemodel = ft_transform_geometry(transform, sourcemodel);
	save(fullfile(dirname, 'ctf', strcat(char(sourcefile), '.mat')), ...
       'sourcemodel');
end
end %transform_to_ctf
