function align_nonlinear(input_img, scalp)
%
% Input:        input_img <string> fullpath to folder, where the aligned
%                                  electrodes shall be stored.
%               scalp <struct> triangular mesh consisting of node positions
%                              and triangle list.
%
% This function aligns the template electrode positions and fiducials from the
% fieldtrip_mni_templates folder to the individual head as defined by the
% input_img and the surface scalp mesh. It does the following for each
% template:
% 1. Loads the template.
% 2. Applies the nonlinear transformation.
% 3. Finetunes the alignment by projecting the points onto the nearast mesh 
%    point of the segmented scalp.
% 4. Saves the aligned electrodes/fiducials.
%
% Example:
% align_nonlinear('/home/user/example_T1.nii', scalp_mesh)
% 
% (c) Nils Harmening, Daniel Miklody, May 2020
% Neurotechnology group, Technische Universität Berlin, Germany

[dirname, base_filename, ext] = fileparts(input_img);
cwd = fileparts(mfilename('fullpath'));

%% SPM nonlinear transformation matrix
nonlin_map = fullfile(dirname, strcat('y_', base_filename, ext));

%% Electrodes
load(fullfile(fileparts(cwd), 'fieldtrip_mni_templates', 'standard_1005.mat'));
elec_aligned = applyDeformationElecs(nonlin_map, elec);
%ft_plot_mesh(scalp); ft_plot_sens(ft_convert_units(elec_aligned, 'mm'));
cfg = [];
cfg.method = 'project';
cfg.headshape = scalp;
elec_aligned = ft_electroderealign(cfg, elec_aligned);
%ft_plot_mesh(scalp); ft_plot_sens(ft_convert_units(elec_aligned, 'mm'));
save(fullfile(dirname, strcat('elec_', 'standard_1005.mat')), 'elec_aligned');

load(fullfile(fileparts(cwd), 'fieldtrip_mni_templates', 'standard_1010.mat'));
elec_aligned = applyDeformationElecs(nonlin_map, elec);
%ft_plot_mesh(scalp_mm); ft_plot_sens(ft_convert_units(elec_aligned, 'mm'));
cfg = [];
cfg.method = 'project';
cfg.headshape = scalp;
elec_aligned = ft_electroderealign(cfg, elec_aligned);
%ft_plot_mesh(scalp_mm); ft_plot_sens(ft_convert_units(elec_aligned, 'mm'));
save(fullfile(dirname, strcat('elec_', 'standard_1010.mat')), 'elec_aligned');

load(fullfile(fileparts(cwd), 'fieldtrip_mni_templates', 'standard_1020.mat'));
elec_aligned = applyDeformationElecs(nonlin_map, elec);
%ft_plot_mesh(scalp_mm); ft_plot_sens(ft_convert_units(elec_aligned, 'mm'));
cfg = [];
cfg.method = 'project';
cfg.headshape = scalp;
elec_aligned = ft_electroderealign(cfg, elec_aligned);
%ft_plot_mesh(scalp_mm); ft_plot_sens(ft_convert_units(elec_aligned, 'mm'));
save(fullfile(dirname, strcat('elec_', 'standard_1020.mat')), 'elec_aligned');

%% Fiducials
load(fullfile(fileparts(cwd), 'fieldtrip_mni_templates', 'fiducials.mat'));
fiducials_aligned = applyDeformationElecs(nonlin_map, fiducials);
%ft_plot_mesh(scalp_mm); ft_plot_sens(ft_convert_units(fiducials_aligned, 'mm'));
cfg = [];
cfg.method = 'project';
cfg.headshape = scalp;
fiducials_aligned = ft_electroderealign(cfg, fiducials_aligned);
%ft_plot_mesh(scalp_mm); ft_plot_sens(ft_convert_units(fiducials_aligned, 'mm'));
save(fullfile(dirname, 'fiducials_nonlinear.mat'), 'fiducials_aligned');

end %align_nonlinear




function [new_elec] = applyDeformationElecs(deffilename, elec)
  % Applies the deformation to electrode positions
  orig_unit=elec.unit;
  new_elec=ft_convert_units(elec,'mm');

  new_pos = applyDeformation(deffilename, new_elec.elecpos);

  new_elec.elecpos=new_pos;
  new_elec.chanpos=new_pos;
  new_elec=ft_convert_units(new_elec,orig_unit);
end %applyDeformationElecs


function [new_sourcemodel] = applyDeformationSourcemodel(deffilename, ...
                                                         sourcemodel)
  % Applies the deformation to sourcemodel positions
  orig_unit = sourcemodel.unit;
  new_sourcemodel=ft_convert_units(sourcemodel, 'mm');

  new_pos = applyDeformation(deffilename, new_sourcemodel.pos);

  new_soucemodel.pos = new_pos;
  new_sourcemodel = ft_convert_units(new_sourcemodel, orig_unit);
end %applyDeformationSourcemodel


function [new_grid] = applyDeformationGrid(deffilename, grid)
  % Applies the deformation to grids
  orig_unit=grid.unit;
  new_grid=ft_convert_units(grid, 'mm');

  new_pos = applyDeformation(deffilename, grid.pos);

  new_grid.pos=new_pos;
  new_grid=ft_convert_units(new_grid,orig_unit);
  if isfield(grid,'norms')
        new_grid.norms=om_normals(new_grid.pos,new_grid.tri);
  end
end %applyDeformationGrid
