function [mesh] = prepare_mesh_iso2mesh(seg, tissue, maxvoxelvol, dofix, ...
                                        downsample)
% 
% Input:        seg <nx x ny x nz matrix> segmented MRI scan
%               tissue <list> of tissue types, e.g. {'brain', 'skull', 'scalp'},
%                             or something like [1 2 3]
%               maxvoxelvol <float> target maximum tetrahedral elem volume
%                                   (default=1)
%               dofix <logical> 0: no mesh repair, 1: mesh repair
%               downsample <integer> downsampling of input seg data (default=1,
%                                    i.e. no downsampling), see fieldtrip
%                                    function ft_volumedownsample
% Output;       mesh <struct> of volumetric tetrahedra mesh
%
% This function prepares the call of vol2mesh for the input of multi-tissue, 
% segmented MRI data. It is based on ft_prepare_mesh from fieldtrip toolbox.
% It does the following:
% 1. Converts and combines seg data into an indexed matrix, where each index is
%    representing one tissue type.
% 2. Calls cgalmesh routine of vol2mesh
% 3. Collects the output as returning mesh struct.
%
% Example:
% mesh = prepare_mesh_iso2mesh(segmentedmri, {'air', scalp','skull','csf',...
%                                             'gray','white'}, 1); 
% 
% (c) Nils Harmening, May 2020
% Neurotechnology group, Technische Universität Berlin, Germany

if nargin > 2
  maxvol = maxvoxelvol;
else
  maxvol = 1;       % default resolution
end
if nargin < 4
  dofix = 1;        %  defult 1: mesh repair
end
if nargin < 5;
  downsample = 1;   % default is no downsampling
end

%% optionally downsample tissue segmentation with fieldtrip
if downsample~=1
  tmpcfg = []
  tmpcfg.downsample = downsample;
  seg = ft_volumedownsample(tmpcfg, seg);
end

%% if tissue not provided, try to get segmented tissue names from segmentedmri
if isempty(tissue)
  seg = ft_datatype_segmentation(seg, 'segmentationstyle', 'indexed');
  fn = fieldnames(seg);
  for i = 1:numel(fn)
    if (numel(seg.(fn{i})) == prod(seg.dim)) && (~strcmp(fn{i}, 'inside'))
      segfield = fn{i};
    end
  end
  tissue = setdiff(unique(seg.(segfield)(:)), 0);
end 
if isempty(tissue)
  error('no tissue provided')
end

if ischar(tissue)
  tissue = {tissue};
end

if iscell(tissue)
  % the code below assumes that it is a probabilistic representation
  if any(strcmp(tissue, 'brain'))
    seg = ft_datatype_segmentation(seg, 'segmentationstyle', ...
                                   'probabilistic', 'hasbrain', 'yes');
  else
    seg = ft_datatype_segmentation(seg, 'segmentationstyle', 'probabilistic');
  end
  % combine all tissue types
  seg_ints = uint8(false(seg.dim));
  for i=1:numel(tissue)
    seg_ints = seg_ints + uint8(seg.(tissue{i}) * i);
    % Check for overlapping / double assinged tissue in segmentation
    % (should not be possible, see preprocessing steps)
    if seg_ints > i
      error('overlapping tissues in input data')
    end
  end
else
  % the code below assumes that it is an indexed representation
  seg = ft_datatype_segmentation(seg, 'segmentationstyle', 'indexed');
  seg_ints = uint8(false(seg.dim));
  % combine all tissue types
  seg_ints = uint8(seg);
end

% check for exisiting iso2mesh toolbox
ft_hastoolbox('iso2mesh', 1);

opt = 2; %2 for FEM % 5 for fNIRS
isovalues = numel(tissue);
[node, elem, face] = vol2mesh(seg_ints, 1:seg.dim(1), 1:seg.dim(2), ...
                              1:seg.dim(3), opt, maxvol, dofix, 'cgalmesh');

[no, fa] = removeisolatednode(node,face(:,1:3));
face = [fa face(:,4)];
[no, el] = removeisolatednode(node,elem(:,1:4));
newelem = meshreorient(no(:,1:3),el(:,1:4));
node = no;
elem = [newelem elem(:,5)];

mesh = [];
mesh.pos = ft_warp_apply(seg.transform, node(:,1:3));
mesh.poslabel = node(:,4);
%mesh.tet = elem(:,[1 2 4 3]); % re-order elements % NO
mesh.tet = elem(:,1:4);
mesh.tetlabel = elem(:,5);
mesh.tri = face(:,1:3);
mesh.trilabel = face(:,4);
mesh.tissues = tissue;

% copy the geometrical units from the input to the output
if isfield(seg, 'unit')
  mesh.unit = seg.unit;
end

% copy the coordinate system from the input to the output
if isfield(seg, 'coordsys')
  mesh.coordsys = seg.coordsys;
end

end %prepare_mesh_iso2mesh
