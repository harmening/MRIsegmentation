function bnd = prepare_nested_meshes(cfg, seg)
% PREPARE_NESTED_MESHES  Boundary surfaces of nested, solid compartments.
%
% ft_prepare_mesh converts a segmentation with several tissues into an indexed
% representation, where every voxel belongs to exactly one tissue. For the
% nested compartments this pipeline builds (scalp contains skull contains csf
% contains cortex) that turns each compartment into a thin layer, and the
% resulting "surface" is the boundary of that layer - a double-walled shell -
% unless FieldTrip's volumefillholes manages to fill the interior again. It
% fails whenever the cavity leaks, e.g. where the skull reaches the edge of
% the field of view (Colin27), and the defect is silent: the mesh is closed
% and consistently oriented, it just encloses the layer instead of the head.
%
% Meshing one tissue at a time keeps each compartment solid, so every surface
% bounds the whole compartment.
%
% Input:   cfg <struct> as for ft_prepare_mesh (tissue, numvertices, method, ...)
%          seg <struct> segmentation with one field per tissue
%
% (c) Eike Middell, October 2026

tissues = cfg.tissue;
numvertices = cfg.numvertices;
if numel(numvertices) == 1
  numvertices = repmat(numvertices, 1, numel(tissues));
end

for i = 1:numel(tissues)
  single = struct('dim', seg.dim, 'transform', seg.transform, 'unit', seg.unit);
  if isfield(seg, 'coordsys')
    single.coordsys = seg.coordsys;
  end
  single.(tissues{i}) = seg.(tissues{i});

  cfgi = cfg;
  cfgi.tissue = tissues(i);
  cfgi.numvertices = numvertices(i);
  b = ft_prepare_mesh(cfgi, single);
  if i == 1
    bnd = b;
  else
    bnd(i) = b; %#ok<AGROW>
  end
end
end %prepare_nested_meshes
