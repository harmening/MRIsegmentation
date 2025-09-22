function create_surface_meshes_cedalion(input_img, numvertices, input_coordsys, ...
                                        output_coordsys)
%
% Input:        input_img <string> fullpath to folder, where the results of
%                                  mysegment.m are stored
%               numvertices <integer> number of desired mesh vertices.
%               input_coordsys <string> name of coordinate system of input_img
%                                       (optional).
%               output_coordsys <string> name of desired coordinate system of
%                                        created meshes (optional).
%
% This function creates triangular surfaces meshes from scalp, skull, csf and
% white and gray matter (=whitegray) on the output of mysegment.m. It does the
% following:
% 1. Loads the tissue masks and combines gray and white matter to whitegray.
% 2. Saves the final segmentation as segmentedmri.
% 3. Calls the cgalsurf algorithm for creating a surface mesh for each tissue
%    class.
% 4. Corrects for segmentation errors and smoothes the surfaces meshes.
% 5. Saves the final meshes.
%
% Example:
% create_surface_meshes('/tmp/head/1/'); expects outputs in directory
% /tmp/head/1/ with input files named: mask_air.nii, mask_skin.nii,
% mask_bone.nii, mask_csf.nii, mask_gray.nii, mask_white.nii
%
% (c) Nils Harmening, Daniel Miklody, May 2020
% Neurotechnology group, Technische Universität Berlin, Germany

[dirname, base_filename, ext] = fileparts(input_img);
if nargin > 1
  numverts = numvertices;
else
  numverts = 1922;
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
if ~numel(dir(fullfile(dirname,'segmentedmri.mat')))
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
    segmentedmri.whitegray = double(segmentedmri.gray|segmentedmri.white);
    segmentedmri.brain = segmentedmri.whitegray|segmentedmri.csf;
    se = ones(3,3,3);
    segmentedmri.csf = segmentedmri.csf|imdilate(double(segmentedmri.brain),se);
    segmentedmri.skull = double((imdilate(double(segmentedmri.csf),se))|...
                                segmentedmri.skull&~segmentedmri.csf);
    segmentedmri.scalp = segmentedmri.scalp|...
                         imdilate(double(segmentedmri.skull),se);
    segmentedmri.skull(:,:,1) = 0;
    segmentedmri.csf(:,:,1:2) = 0;
    segmentedmri.whitegray(:,:,1:3) = 0;

    fields = {'air','brain'};
    segmentedmri = rmfield(segmentedmri,fields);
    save(fullfile(dirname,'segmentedmri'), 'segmentedmri');
end
clear c1 c2 c3 c4 c5 c6

if nargin > 3
  segmentedmri.anatomy = mri.anatomy;
  segmentedmri = ft_convert_coordsys(segmentedmri, output_coordsys);
  save(fullfile(dirname,strcat('segmentedmri_',output_coordsys)), ...
       'segmentedmri');
end

%% Build surface meshes (4 meshes)
if ~numel(dir(fullfile(dirname, strcat('bnd4_', num2str(numverts), '.mat'))))
    load(fullfile(dirname,'segmentedmri'));
    fields = {'white', 'gray'};
    segmentedmri = rmfield(segmentedmri,fields);
    cfg2 = [];
    cfg2.tissue = {'scalp','skull','csf','whitegray'};
    cfg2.numvertices = [numverts numverts numverts numverts numverts];
    cfg2.transform = segmentedmri.transform;
    cfg2.method = 'iso2mesh'; %'projectmesh' 

    bnd = ft_prepare_mesh(cfg2, segmentedmri);
    save(fullfile(dirname, strcat('bnd4_', num2str(numverts))), 'bnd');
    % cedalion export
    scalp.pos = bnd(1).pos;
    scalp.tri = bnd(1).tri;
    invtransform = inv(mri.transform);
    scalp = ft_transform_geometry(invtransform, scalp);
    om_save_stl(fullfile(dirname, strcat('mask_scalp', num2str(numverts), ...
                                         '.stl')), scalp.pos, scalp.tri);

    whitegray.pos = bnd(4).pos;
    whitegray.tri = bnd(4).tri;
    invtransform = inv(mri.transform);
    whitegray = ft_transform_geometry(invtransform, whitegray);
    om_save_stl(fullfile(dirname, strcat('mask_brain', num2str(numverts), ...
                                         '.stl')), whitegray.pos, whitegray.tri);

end
clear segmentedmri

%% Correct outliers and smooth skull&brain&csf
if ~numel(dir(fullfile(dirname, strcat('bnd4_', num2str(numverts), ...
                       '_corrected.mat'))))
    load(fullfile(dirname, strcat('bnd4_', num2str(numverts))));
    error_threshold = 25;
    smooth_meshes = [1 2 3 4];
    new_bnd = correct_bnd_errors(bnd, error_threshold, smooth_meshes);
    new_bnd = ft_convert_units(new_bnd,'m');
    save(fullfile(dirname, strcat('bnd4_', num2str(numverts), ...
                                  '_corrected')), 'new_bnd');
    om_save_tri(fullfile(dirname, strcat('bnd4_', num2str(numverts), ...
                                         '_corrected_scalp.tri')), ...
                new_bnd(1).pos, new_bnd(1).tri);
    %om_save_tri(fullfile(dirname, strcat('bnd4_', num2str(numverts), ...
    %                                     '_corrected_skull.tri')), ...
    %            new_bnd(2).pos, new_bnd(2).tri);
    %om_save_tri(fullfile(dirname, strcat('bnd4_', num2str(numverts), ...
    %                                     '_corrected_csf.tri')), ...
    %            new_bnd(3).pos, new_bnd(3).tri);
    %om_save_tri(fullfile(dirname, strcat('bnd4_', num2str(numverts), ...
    %                                     '_corrected_cortex.tri')), ...
    %            new_bnd(4).pos, new_bnd(4).tri);
    % cedalion export
    new_bnd = ft_convert_units(new_bnd,'mm');
    scalp.pos = new_bnd(1).pos;
    scalp.tri = new_bnd(1).tri;
    invtransform = inv(mri.transform);
    scalp = ft_transform_geometry(invtransform, scalp);
    om_save_stl(fullfile(dirname, strcat('mask_scalp', num2str(numverts), ...
                                         '_smoothed.stl')), ...
                scalp.pos, scalp.tri);
    write_obj(fullfile(dirname, 'mask_scalp.obj'), scalp.pos, scalp.tri);

end

end %create_surface_meshes




function [nrm] = normals(pnt, tri, opt)
  % Compute the surface normals of a triangular mesh for each triangle or for
  % each vertex.
  % Adapted from fieldtrip/private/normals.m

	npnt = size(pnt,1);
	ntri = size(tri,1);

	% shift to center
	pnt(:,1) = pnt(:,1)-mean(pnt(:,1),1);
	pnt(:,2) = pnt(:,2)-mean(pnt(:,2),1);
	pnt(:,3) = pnt(:,3)-mean(pnt(:,3),1);

	% compute triangle normals
	v2 = pnt(tri(:,2),:) - pnt(tri(:,1),:);
	v3 = pnt(tri(:,3),:) - pnt(tri(:,1),:);
	nrm_tri = cross(v2, v3);
	if strcmp(opt, 'faces')
	  % normalise
	  nrm = nrm_tri ./ (sqrt(sum(nrm_tri.^2, 2)) * ones(1,3));

        elseif strcmp(opt, 'vertices')
	  % compute vertex normals
	  nrm_pnt = zeros(npnt, 3);
	  for i=1:ntri
	    nrm_pnt(tri(i,1),:) = nrm_pnt(tri(i,1),:) + nrm_tri(i,:);
	    nrm_pnt(tri(i,2),:) = nrm_pnt(tri(i,2),:) + nrm_tri(i,:);
	    nrm_pnt(tri(i,3),:) = nrm_pnt(tri(i,3),:) + nrm_tri(i,:);
	  end
	  % normalise
	  nrm = nrm_pnt ./ (sqrt(sum(nrm_pnt.^2, 2)) * ones(1,3));
	end
end %normals


function [c] = cross(a,b)
  % Fast cross product to replace the MATLAB standard version
  c = [a(:,2).*b(:,3)-a(:,3).*b(:,2) a(:,3).*b(:,1)-a(:,1).*b(:,3) ...
       a(:,1).*b(:,2)-a(:,2).*b(:,1)];
end %cross


function [bnd] = correct_bnd_errors(bnd, errorthreshold, smooth)
  % Correction of triangular surface meshes for abnormal vertex outliers and
  % application of smoothing routines.
  num_corrected=0;
  for ii=1:numel(bnd)
      for j=1:size(bnd(ii).tri,1)
          if (norm(bnd(ii).pos(bnd(ii).tri(j,1),:) - ...
                   bnd(ii).pos(bnd(ii).tri(j,2),:)) > errorthreshold)
              if (norm(bnd(ii).pos(bnd(ii).tri(j,1),:) - ...
                       bnd(ii).pos(bnd(ii).tri(j,3),:)) > errorthreshold)
                  [nb_row nb_col] = find(bnd(ii).tri == bnd(ii).tri(j,1));
                  nbs = unique(bnd(ii).tri(nb_row,:));
                  nbs(nbs==bnd(ii).tri(j,1)) = [];
                  bnd(ii).pos(bnd(ii).tri(j,1),:) = sum(bnd(ii).pos(nbs,:)) ...
                                                    / numel(nbs);
                  num_corrected = num_corrected+1;
              elseif (norm(bnd(ii).pos(bnd(ii).tri(j,2),:) - ...
                           bnd(ii).pos(bnd(ii).tri(j,3),:)) > errorthreshold)
                  [nb_row nb_col] = find(bnd(ii).tri == bnd(ii).tri(j,2));
                  nbs = unique(bnd(ii).tri(nb_row,:));
                  nbs(nbs==bnd(ii).tri(j,2)) = [];
                  bnd(ii).pos(bnd(ii).tri(j,2),:) = sum(bnd(ii).pos(nbs,:)) ...
                                                    / numel(nbs);
                  num_corrected = num_corrected+1;
              end
          elseif (norm(bnd(ii).pos(bnd(ii).tri(j,2),:) - ...
                       bnd(ii).pos(bnd(ii).tri(j,3),:)) > errorthreshold) ...
                 && (norm(bnd(ii).pos(bnd(ii).tri(j,1),:) - ...
                         bnd(ii).pos(bnd(ii).tri(j,3),:)) > errorthreshold)
              [nb_row nb_col] = find(bnd(ii).tri == bnd(ii).tri(j,3));
              nbs = unique(bnd(ii).tri(nb_row,:));
              nbs(nbs==bnd(ii).tri(j,3)) = [];
              bnd(ii).pos(bnd(ii).tri(j,3),:) = sum(bnd(ii).pos(nbs,:)) / ...
                                                numel(nbs);
              num_corrected = num_corrected+1;
          end
      end
  end

  if smooth
      for ii=smooth
          bnd(ii).pos = lpflow_trismooth(bnd(ii).pos,bnd(ii).tri);
      end
  end
  disp(strcat('Corrected ', num2str(num_corrected), ' vertex positions.'));
end %correct_bnd_errors


function xyzn=lpflow_trismooth(xyz,t)
  % Laplace flow mesh smoothing for vertex ring
  %
  % Reference:    1) Zhang and Hamza, (2006) Vertex-based anisotropic
  %               smoothing of 3D mesh data, IEEE CCECE
  % Acknowledgement:
  %               Q. Fang: iso2mesh (http://iso2mesh.sf.net)
  %
  % Input:        xyz <nx3> vertex coordinates
  %               t <mx3> triangulation index array
  % Output:       xyzn <nx3> updates vertex coordinates
  % Version:      1
  % JOK 300709

  % I/O check:
  if nargin~=2
      error('Wrong # of input')
  end
  if nargout ~= 1
      error('Output is a single array, wrong designation!')
  end
  nt= size(t);
  if nt(2)~=3
      error('Triangle element matrix should be mx3!')
  end
  mp = size(xyz);
  if mp(2)~=3
      error('Vertices should be nx3!')
  end

  % Initialize etc.
  k=0;
  [conn,connnum,count]=neighborelem(t,max(max(t)));
  xyzn = xyz;
  while k<length(xyz)
      k=k+1;
      % Find vertex neighborhood
      indt01=conn{k}; % Element indices
      indv01 = t(indt01,:);
      indv01 = unique(indv01(:));
      vdist = xyz(indv01,:)-repmat(xyz(k,:),length(indv01),1);
      dist = sqrt(sum(vdist.*vdist,2));
      indaux1 = find(dist==0);
      vdist(indaux1,:) =[];
      d=length(vdist); % Cardinality of vertex
      if isempty(dist)
          xyzn(k,:)  = NaN;
      end
      vcorr = sum(vdist/d,1);
      xyzn(k,:) = xyz(k,:)+vcorr;
  end
end % lpflow_trismooth


function om_save_stl(filename, pos, tri, solidname)
	% Save a triangular mesh to ASCII STL (facet normals computed from geometry).
	if nargin < 4 || isempty(solidname)
			[~,solidname,~] = fileparts(filename);
	end
	fid = fopen(filename, 'w');
	if fid < 0, error('Cannot open %s for writing.', filename); end
	fprintf(fid, 'solid %s\n', solidname);

	v1 = pos(tri(:,2),:) - pos(tri(:,1),:);
	v2 = pos(tri(:,3),:) - pos(tri(:,1),:);
	n  = cross(v1, v2);                         % uses the vectorized cross below
	nn = sqrt(sum(n.^2,2)); n = n ./ (nn + eps);

	for i = 1:size(tri,1)
			fprintf(fid, '  facet normal %.6g %.6g %.6g\n', n(i,1), n(i,2), n(i,3));
			fprintf(fid, '    outer loop\n');
			fprintf(fid, '      vertex %.6g %.6g %.6g\n', pos(tri(i,1),1), pos(tri(i,1),2), pos(tri(i,1),3));
			fprintf(fid, '      vertex %.6g %.6g %.6g\n', pos(tri(i,2),1), pos(tri(i,2),2), pos(tri(i,2),3));
			fprintf(fid, '      vertex %.6g %.6g %.6g\n', pos(tri(i,3),1), pos(tri(i,3),2), pos(tri(i,3),3));
			fprintf(fid, '    endloop\n');
			fprintf(fid, '  endfacet\n');
	end
	fprintf(fid, 'endsolid %s\n', solidname);
	fclose(fid);
end


function write_obj(path, V, F)
    fid = fopen(path,'w'); assert(fid>0, 'Cannot write OBJ');
    for i=1:size(V,1), fprintf(fid, 'v %.6f %.6f %.6f\n', V(i,1), V(i,2), V(i,3)); end
    for i=1:size(F,1), fprintf(fid, 'f %d %d %d\n', F(i,1), F(i,2), F(i,3)); end
    fclose(fid);
end
