function transform_to_ctf(input_img)
%
% Input:        input_img <string> fullpath to T1 MRI image to be segmented
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


%% Convert meshes
cfg = [];
cfg.method        = 'fiducial';
cfg.coordsys      = 'ctf';
cfg.fiducial.nas  = nas; %position of NAS
cfg.fiducial.lpa  = lpa; %position of LPA
cfg.fiducial.rpa  = rpa; %position of RPA

new_bnd = [];
new_bnd.coordsys = 'ctf';
if numel(dir(fullfile(dirname, 'bnd4_1922_corrected.mat')))
  load(fullfile(dirname, 'bnd4_1922_corrected.mat'));
  shells = {'scalp', 'skull', 'csf', 'cortex'};
  for i=1:4
    new_bnd_i = ft_meshrealign(cfg, new_bnd(i));
    new_bnd(i).pos = new_bnd_i.pos;
    new_bnd(i).tri = new_bnd_i.tri;
    om_save_tri(fullfile(dirname, 'ctf', strcat('bnd4_1922_corrected_', ...
                                                char(shells(i)), '.tri')), ...
                new_bnd(i).pos, new_bnd(i).tri);
  end
end
save(fullfile(dirname, 'ctf', 'bnd4_1922_corrected.mat'), 'new_bnd');

if numel(dir(fullfile(dirname, 'mesh6_maxvoxvol2.mat')))
  load(fullfile(dirname, 'mesh6_maxvoxvol2.mat'));
  transform = ft_headcoordinates(nas, lpa, rpa, 'ctf');
  mesh.pos = ft_transform_geometry(transform, mesh.pos);
  % same as 
  %mesh = ft_meshrealign(cfg, mesh);
  save(fullfile(dirname, 'ctf', strcat('mesh6_maxvoxvol2.mat')), 'mesh');
end

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
sourcefiles = {'sourcemodel4000'};
for sourcefile = sourcefiles
	load(fullfile(dirname, strcat(char(sourcefile), '.mat')));
	sourcemodel = ft_transform_geometry(transform, sourcemodel);
	save(fullfile(dirname, 'ctf', strcat(char(sourcefile), '.mat')), ...
       'sourcemodel');
end
end %transform_to_ctf
