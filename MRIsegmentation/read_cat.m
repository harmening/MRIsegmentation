function read_cat(cat12_path, input_img, num_sources, num_cortex_verts);
%
% Input:        cat12_path <string> fullpath to cat12 installation.
%               input_img <string> fullpath to T1 MRI image to be segmented by
%                                  CAT12.
%               num_sources <integer> number of sourcemodel points.
%
% This function reads the output of a CAT12 segmentation and creates a realistic
% sourcemodel from it. It does the following:
% 1. Loads central, sphere and thickness mesh files for each hemisphere.
% 2. Downsamples the read meshes up to the desired sourcemodel resolution.
% 3. Merge left and right hemisphere meshes.
% 4. Calculates sourcemodel points at 2/3 distance between gm and wm.
% 5. Calculate sourcepoint orientations as surface normals.
%
% Example:
% read_cat('/usr/local/cat12', '/tmp/head/1/T1.nii', 2000);
% Created sourcemodel is stored at the location of the input_img.
%
% (c) Nils Harmening, May 2020
% Neurotechnology group, Technische Universität Berlin, Germany

if nargin < 4 || isempty(num_cortex_verts)
      num_cortex_verts = 60000;
end

[dirname, name, ext] = fileparts(input_img);

cortexmesh_fn_mat  = fullfile(dirname, sprintf('mask_brain%d.mat',  num_cortex_verts));
cortexmesh_fn_tri  = fullfile(dirname, sprintf('mask_brain%d.tri',  num_cortex_verts));
cortexmesh_fn_stl  = fullfile(dirname, sprintf('mask_brain%d.stl',  num_cortex_verts));

if (~numel(dir(fullfile(dirname, strcat('sourcemodel', num2str(num_sources), ...
                                       '.mat')))) || ...
    ~numel(dir(cortexmesh_fn_stl)))

  %% Important CAT12 output files
  fn_tess_lh = fullfile(dirname, 'surf', strcat('lh.central.', name, '.gii'));
  fn_tess_rh = fullfile(dirname, 'surf', strcat('rh.central.', name, '.gii'));
  fn_tess_lsph = fullfile(dirname, 'surf', strcat('lh.sphere.', name, '.gii'));
  fn_tess_rsph = fullfile(dirname, 'surf', strcat('rh.sphere.', name, '.gii'));
  fn_thick_lh = fullfile(dirname, 'surf', strcat('lh.thickness.', name));
  fn_thick_rh = fullfile(dirname, 'surf', strcat('rh.thickness.', name));

  %% Import CAT surfaces
  % left pial
  if ~isempty(fn_tess_lh)
      bnd_lh = import_gii_surfaces(fn_tess_lh);
      num_verts_orig_l = size(bnd_lh.pos, 1);
      save(fullfile(dirname, 'surf', 'lh_surf.mat'), 'bnd_lh');
      om_save_tri(fullfile(dirname, 'surf', 'lh_surf.tri'), bnd_lh.pos, bnd_lh.tri);
  end
  % left corresponding sphere mesh
  if ~isempty(fn_tess_lsph)
      lh_sphere = import_gii_surfaces(fn_tess_lsph);
      save(fullfile(dirname, 'surf', 'lh_sphere.mat'), 'lh_sphere');
      om_save_tri(fullfile(dirname, 'surf', 'lh_sphere.tri'), lh_sphere.pos, lh_sphere.tri);
  end
  % right pial
  if ~isempty(fn_tess_rh)
      bnd_rh = import_gii_surfaces(fn_tess_rh);
      num_verts_orig_r = size(bnd_rh.pos, 1);
      save(fullfile(dirname, 'surf', 'rh_surf.mat'), 'bnd_rh');
      om_save_tri(fullfile(dirname, 'surf', 'rh_surf.tri'), bnd_rh.pos, bnd_rh.tri);
  end
  % right corresponding sphere mesh
  if ~isempty(fn_tess_rsph)
      rh_sphere = import_gii_surfaces(fn_tess_rsph);
      save(fullfile(dirname, 'surf', 'rh_sphere.mat'), 'rh_sphere');
      om_save_tri(fullfile(dirname, 'surf', 'rh_sphere.tri'), rh_sphere.pos, rh_sphere.tri);
  end

  %% cedalion export (detailed CAT12 cortical surface mesh (NOT shifted inward))
  if exist('bnd_lh','var') && exist('bnd_rh','var') && ...
     (~exist(cortexmesh_fn_mat,'file') || ~exist(cortexmesh_fn_tri,'file'))

      num_vertshemi_cortex = max(1, round(num_cortex_verts/2));
      disp(['CORTEX MESH> Downsampling left hemi to ~', num2str(num_vertshemi_cortex), ' vertices.']);
      [cort_lh_low, ~, ~] = tess_downsize(bnd_lh, num_vertshemi_cortex);

      disp(['CORTEX MESH> Downsampling right hemi to ~', num2str(num_vertshemi_cortex), ' vertices.']);
      [cort_rh_low, ~, ~] = tess_downsize(bnd_rh, num_vertshemi_cortex);

      cortex_mesh = tess_concatenate({cort_lh_low, cort_rh_low}, ...
                      sprintf('cortex_%dV', size(cort_lh_low.pos,1) + size(cort_rh_low.pos,1)), 'Cortex');

      % Save as one merged, clean mesh
      %save(cortexmesh_fn_mat, 'cortex_mesh');
      %om_save_tri(cortexmesh_fn_tri, cortex_mesh.pos, cortex_mesh.tri);
      mri = ft_read_mri(input_img);
      invtransform = inv(mri.transform);
      cortex_mesh = ft_transform_geometry(invtransform, cortex_mesh);
      om_save_stl(cortexmesh_fn_stl, cortex_mesh.pos, cortex_mesh.tri);
  end





  %% Import CAT12 cortical thicknesses
  if ~isempty(fn_thick_lh) && ~isempty(fn_thick_lh)
      % right pial
      thick_rh = read_curv(fn_thick_rh);
      pos = bnd_rh.pos;
      tri = bnd_rh.tri;
      nrms = vert_normals(pos, tri);
      for i=1:length(pos)
        pos(i,:) = pos(i,:) - nrms(i,:) * thick_rh(i)*1/3;
      end
      bnd_rh.pos = pos;
      % left pial
      thick_lh = read_curv(fn_thick_lh);
      pos = bnd_lh.pos;
      tri = bnd_lh.tri;
      nrms = vert_normals(pos, tri);
      for i=1:length(pos)
        pos(i,:) = pos(i,:) - nrms(i,:) * thick_lh(i)*1/3;
      end
      bnd_lh.pos = pos;
  end

  %% Downsample left and right hemispheres
  num_vertshemi = round(num_sources / 2);
  if ~isempty(fn_tess_rh)
      disp(strcat('Import CAT12 folder', 'Downsampling right pial.'));
      [tess_rh_low, i_rh_low, x_rh_low] = tess_downsize(bnd_rh, num_vertshemi);
  end
  if ~isempty(fn_tess_lh)
      disp(strcat('Import CAT12 folder', 'Downsampling left pial.'));
      [tess_lh_low, i_lh_low, x_lh_low] = tess_downsize(bnd_lh, num_vertshemi);
  end


  %% Merge hemispheres
  if ~isempty(fn_tess_lh) && ~isempty(fn_tess_rh)
      % Hi-resolution surface (not needed for now):
      cortex_high  = tess_concatenate({bnd_lh, bnd_rh}, ...
                                   sprintf('cortex_%dV', num_verts_orig_l ...
                                           + num_verts_orig_r), 'Cortex');
      cortex_low = tess_concatenate({tess_lh_low, tess_rh_low}, ...
                                   sprintf('cortex_%dV', length(x_lh_low) + ...
                                                         length(x_rh_low)), ...
                                           'Cortex');
  else
      coerex_high_file = [];
      cortex_low_file = [];
  end

  %% Calculate surfaces
  pos = cortex_low.pos;
  tri = cortex_low.tri;
  nrms = vert_normals(pos, tri);
  if (sum(sum(isnan(nrms))) > 0)
    % Correct normals by using the average normal of neighboring vertices
    [row, col] = find(isnan(nrms));
    nan_idx = unique(row);
    for idx = 1:length(nan_idx)
      [row, col] = find(tri == nan_idx(idx));
      neighbors = [];
      t = unique(row);
      if sum(t) == 0
        disp('Found no tris!');
        disp(t);
      end
    for t_idx = 1:length(t)
	    neighbors = [neighbors; tri(t(t_idx),1); tri(t(t_idx),2); tri(t(t_idx),3)];
    end
    neighbors = unique(neighbors(neighbors~=nan_idx(idx)));
    nrm = nansum(nrms(neighbors,:),1);
    nrm = nrm / norm(nrm);
    nrms(nan_idx(idx),:) = nrm;
    end
  end
  sourcemodel = [];
  sourcemodel.pos = pos;
  sourcemodel.tri = tri;
  sourcemodel.vec = nrms;
  sourcemodel.unit = 'mm';
  sourcemodel.inside = true(length(sourcemodel.pos), 1);
  save(fullfile(dirname, strcat('sourcemodel', num2str(num_sources), ...
                                '.mat')), 'sourcemodel');
  om_save_tri(fullfile(dirname, strcat('sourcemodel', num2str(num_sources), ...
                                       '.tri')), sourcemodel.pos, ...
              sourcemodel.tri);


  nrms = sourcemodel.vec;
  if (sum(sum(isnan(nrms))) > 0)
    % Correct normals by using the average normal of neighboring vertices
    [row, col] = find(isnan(nrms));
    nan_idx = unique(row);
    for idx = 1:length(nan_idx)
      [row, col] = find(tri == nan_idx(idx));
      neighbors = [];
      t = unique(row);
      if sum(t) == 0
        disp('found no tris!')
        disp(t)
      end
      for t_idx = 1:length(t)
        neighbors = [neighbors; tri(t(t_idx),1); tri(t(t_idx),2); ...
                     tri(t(t_idx),3)];
      end
      neighbors = unique(neighbors(neighbors~=nan_idx(idx)));
      nrm = nansum(nrms(neighbors,:),1);
      nrm = nrm / norm(nrm);
      nrms(nan_idx(idx),:) = nrm;
    end
    if (sum(sum(isnan(nrms))) > 0)
      disp(strcat('Could not compute ', num2str(sum(sum(isnan(nrms)))), ...
                  'sourcemodel normals in ', dirname))
    end
  end
  sourcemodel.vec = nrms;
  save(fullfile(dirname, strcat('sourcemodel', num2str(num_sources), '.mat')),...
         'sourcemodel');
end
end %read_cat




function new_bnd = import_gii_surfaces(filename)
  % Imports gii surface files
  % Adapted from brainstorm
  g = gifti(filename);
  Tess.pos = double(g.vertices);
  Tess.tri = g.faces(:,[2 1 3]);
  Tess.mat = g.mat;
	% Only one surface
	if (length(Tess) == 1)
			% Surface mesh
			if isfield(Tess, 'tri') % Volume meshes do not have triangulation field
					new_bnd.pos = Tess(1).pos;
					new_bnd.tri = Tess(1).tri;
			% Volume FEM mesh
			else
					new_bnd = Tess;
			end
      [tmp__, fBase, fExt] = fileparts(filename);
      importedBaseName = [fBase, strrep(fExt, '.', '_')];
      importedBaseName = strrep(importedBaseName, 'tess_', '');
      importedBaseName = strrep(importedBaseName, '_tess', '');
      new_bnd.Comment = importedBaseName;
	% Multiple surfaces
	else
		new_bnd = tess_concatenate(Tess, [], []);
		new_bnd.iAtlas  = find(strcmpi({new_bnd.Atlas.Name}, 'Structures'));
	end
end %import_gii_surfaces


function [NewTessMat, I, J] = tess_downsize(bnd, newNbVertices)
  % Downsize the new tesselation
  % Adapted from brainstorm
  NewTessMat = [];
  NewTessMat.Comment = bnd.Comment;
  I = [];
  J = [];
  oldNbVertices = size(bnd.pos, 1);
  if (newNbVertices >= oldNbVertices)
        disp(sprintf('TESS> Surface has ', num2str(oldNbVertices), ...
                     'vertices out of ', num2str(newNbVertices)));
        newNbVertices = oldNbVertices;
  end

	% Prepare variables
	bnd.tri    = double(bnd.tri);
	bnd.pos = double(bnd.pos);
	dsFactor = newNbVertices / size(bnd.pos, 1);
  % Matlab's reducepatch:
  % Reduce number of vertices
  [NewTessMat.tri, NewTessMat.pos] = reducepatch(bnd.tri, bnd.pos, ...
                                                 dsFactor);
  % Find the vertices that were kept by reducepatch
  [tmp, I, J] = intersect(bnd.pos, NewTessMat.pos, 'rows');
  % Re-order the vertices so that they are in the same order in the output surface
  [I, iSort] = sort(I);
  NewTessMat.pos = bnd.pos(I,:);
  J = J(iSort);
  % Re-order the vertices in the faces
  iSortFaces(J) = 1:length(J);
  NewTessMat.tri = iSortFaces(NewTessMat.tri);
  % Set the
  J = (1:length(J))';
end %tess_downsize


function new_bnd = tess_concatenate(bnds, NewComment, fileType)
	% Concatenate tesselletions
  % Adapted from brainstorm
	new_bnd = [];
    new_bnd.pos = [];
    new_bnd.tri = [];
	isLeft = 0;
	isRight = 0;
	isWhite = 0;
	isCortex = 0;
	isAseg = 0;
	isSave = 1;
	for iFile = 1:length(bnds)
    % Load tesselation
    if iscell(bnds)
        old_bnd = bnds{iFile};
        if isempty(old_bnd)
            continue
        end
    % Files already loaded in calling function
    else
        old_bnd = bnds(iFile);
        isSave = 0;
    end

	% Detect if right/left hemisphere
    if ~isempty(strfind(old_bnd.Comment, 'lh.')) || ...
          ~isempty(strfind(old_bnd.Comment, 'Lhemi')) || ...
          ~isempty(strfind(old_bnd.Comment, 'Lwhite')) || ...
          ~isempty(strfind(old_bnd.Comment, 'left'))
        isLeft = 1;
        scoutTag = ' L';
        scoutHemi = 'L';
        scoutComment = 'Cortex';
    elseif ~isempty(strfind(old_bnd.Comment, 'rh.')) || ...
           ~isempty(strfind(old_bnd.Comment, 'Rhemi')) || ...
           ~isempty(strfind(old_bnd.Comment, 'Rwhite')) || ...
           ~isempty(strfind(old_bnd.Comment, 'right'))
        isRight = 1;
        scoutTag = ' R';
        scoutHemi = 'R';
        scoutComment = 'Cortex';
    % Detect based on comment (tag ' L' or ' R' already present)
    elseif (length(old_bnd.Comment) > 2) && ...
           strcmpi(old_bnd.Comment(end-1:end), ' L')
        scoutTag = ' L';
        scoutHemi = 'L';
        scoutComment = old_bnd.Comment;
    elseif (length(old_bnd.Comment) > 2) && ...
            strcmpi(old_bnd.Comment(end-1:end), ' R')
        scoutTag = ' R';
        scoutHemi = 'R';
        scoutComment = old_bnd.Comment;
    % Guess based on the coordinates
    else
        if (nnz(old_bnd.pos(:,2) > 0) > 5 * nnz(old_bnd.pos(:,2) < 0))
            scoutTag = ' L';
            scoutHemi = 'L';
        elseif (nnz(old_bnd.pos(:,2) < 0) > 5 * nnz(old_bnd.pos(:,2) > 0))
            scoutTag = ' R';
            scoutHemi = 'R';
        else
            scoutTag = '';
            scoutHemi = 'U';
        end
        scoutComment = old_bnd.Comment;
    end
		% Detect some specific types of surfaces
    if ~isempty(strfind(old_bnd.Comment, 'white'))
        isWhite = 1;
    end
    if ~isempty(strfind(old_bnd.Comment, 'cortex_'))
        isCortex = 1;
    end
    if ~isempty(strfind(old_bnd.Comment, 'aseg'))
        isAseg = 1;
    end
    % Concatenate current sub-tess to final tesselation structure
    offsetVertices   = size(new_bnd.pos,1);
    new_bnd.tri    = [new_bnd.tri; old_bnd.tri + offsetVertices];
    new_bnd.pos = [new_bnd.pos; old_bnd.pos];



		% Concatenate FreeSurfer registration spheres
    if isfield(old_bnd, 'Reg') && isfield(old_bnd.Reg, 'Sphere') && ...
          isfield(old_bnd.Reg.Sphere, 'pos') && ~isempty(old_bnd.Reg.Sphere.pos)
        if ~isfield(new_bnd, 'Reg') || ~isfield(new_bnd.Reg, 'Sphere') || ...
              ~isfield(new_bnd.Reg.Sphere, 'pos')
            new_bnd.Reg.Sphere.pos = old_bnd.Reg.Sphere.pos;
        else
            new_bnd.Reg.Sphere.pos = [new_bnd.Reg.Sphere.pos; ...
                                      old_bnd.Reg.Sphere.pos];
        end
    end
    % Concatenate BrainSuite registration squares
    if isfield(old_bnd, 'Reg') && isfield(old_bnd.Reg, 'Square') && ...
          isfield(old_bnd.Reg.Square, 'pos') && ~isempty(old_bnd.Reg.Square.pos)
        if ~isfield(new_bnd, 'Reg') || ~isfield(new_bnd.Reg, 'Square') || ...
              ~isfield(new_bnd.Reg.Square, 'pos')
            new_bnd.Reg.Square.pos      = old_bnd.Reg.Square.pos;
            new_bnd.Reg.AtlasSquare.pos = old_bnd.Reg.AtlasSquare.pos;
        else
            new_bnd.Reg.Square.pos      = [new_bnd.Reg.Square.pos; ...
                                           old_bnd.Reg.Square.pos];
            new_bnd.Reg.AtlasSquare.pos = [new_bnd.Reg.AtlasSquare.pos; ...
                                           old_bnd.Reg.AtlasSquare.pos];
        end
    end
	end

	% Detect surfaces types
	if isLeft && isRight
			% File type: Cortex
			if isempty(fileType)
					fileType = 'Cortex';
			end
			% White matter
			if isWhite
					fileTag = 'cortex_white';
					if isempty(NewComment)
							NewComment = sprintf('white_%dV', length(new_bnd.pos));
					end
			% Pial/cortex external envelope
			else
					fileTag = 'cortex_pial';
					if isempty(NewComment)
							NewComment = sprintf('cortex_%dV', length(new_bnd.pos));
					end
			end
	elseif isCortex && isAseg
			fileTag = 'cortex_mixed';
			if isempty(fileType)
					fileType = 'Cortex';
			end
			if isempty(NewComment)
					NewComment = sprintf('cortex_mixed_%dV', length(new_bnd.pos));
			end
	else
			fileTag = 'concat';
			if isempty(fileType)
					fileType = 'Other';
			end
			if isempty(NewComment)
					NewComment = 'New surface';
			end
	end
	% Surface comments
	new_bnd.Comment = NewComment;
end %tess_concatenate




function [nrm] = vert_normals(pnt, tri)
	% Surface normals of a triangular mesh for each triangle or for each vertex
  % Adapted from fieldtrip
	npnt = size(pnt,1);
	ntri = size(tri,1);

	% shift to center
	pnt(:,1) = pnt(:,1)-mean(pnt(:,1),1);
	pnt(:,2) = pnt(:,2)-mean(pnt(:,2),1);
	pnt(:,3) = pnt(:,3)-mean(pnt(:,3),1);

	% compute normal of triangular
	v2 = pnt(tri(:,2),:) - pnt(tri(:,1),:);
	v3 = pnt(tri(:,3),:) - pnt(tri(:,1),:);
	nrm_tri = cross(v2, v3);

  % compute vertex normals
  nrm_pnt = zeros(npnt, 3);
  for i=1:ntri
    nrm_pnt(tri(i,1),:) = nrm_pnt(tri(i,1),:) + nrm_tri(i,:);
    nrm_pnt(tri(i,2),:) = nrm_pnt(tri(i,2),:) + nrm_tri(i,:);
    nrm_pnt(tri(i,3),:) = nrm_pnt(tri(i,3),:) + nrm_tri(i,:);
  end
  % normalise the direction vectors to have length one
  nrm = nrm_pnt ./ (sqrt(sum(nrm_pnt.^2, 2)) * ones(1,3));
end %normals


function [c] = cross(a,b)
	c = [a(:,2).*b(:,3)-a(:,3).*b(:,2) a(:,3).*b(:,1)-a(:,1).*b(:,3) ...
       a(:,1).*b(:,2)-a(:,2).*b(:,1)];
end %cross


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
