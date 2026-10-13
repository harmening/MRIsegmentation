function bnd = close_mesh_holes(bnd, maxloop)
% CLOSE_MESH_HOLES  Close small holes in triangulated surface meshes.
%
% cgalsurf (through ft_prepare_mesh) occasionally returns surfaces with a few
% small holes, in particular when it hits the vertex budget before satisfying
% its angular bound ("reach maximum node limit"). An open surface silently
% breaks everything that needs an inside/outside decision: BEM, point-in-head
% tests and photon simulations with sources placed on the surface.
%
% Each hole is closed with a triangle fan around the centroid of its boundary
% loop, oriented to match the triangles it is attached to. Only the hole is
% touched; the rest of the mesh keeps its vertices and triangles.
%
% Input:   bnd     <struct array> meshes with .pos and .tri
%          maxloop <integer> largest boundary loop to close (default 200
%                            vertices); larger openings are left alone and
%                            reported, as they indicate a real problem rather
%                            than a meshing glitch
%
% (c) Eike Middell, October 2026

if nargin < 2 || isempty(maxloop)
  maxloop = 200;
end

for ii = 1:numel(bnd)
  pos = bnd(ii).pos;
  tri = bnd(ii).tri;

  % directed edges; in a closed, consistently oriented mesh every edge occurs
  % once in each direction, so a boundary edge is one without its reverse
  edges = [tri(:, [1 2]); tri(:, [2 3]); tri(:, [3 1])];
  [~, loc] = ismember(edges(:, [2 1]), edges, 'rows');
  boundary = edges(loc == 0, :);
  if isempty(boundary)
    continue
  end

  % walk the boundary edges into loops: each vertex has one outgoing edge
  next = containers.Map(boundary(:, 1), boundary(:, 2));
  unused = containers.Map(boundary(:, 1), true(size(boundary, 1), 1));
  nclosed = 0;
  nleft = 0;
  while ~isempty(keys(unused))
    k = keys(unused);
    start = k{1};
    loop = start;
    v = next(start);
    remove(unused, start);
    while v ~= start && isKey(unused, v)
      loop(end + 1) = v; %#ok<AGROW>
      nxt = next(v);
      remove(unused, v);
      v = nxt;
    end
    if v ~= start || numel(loop) < 3
      % not a simple cycle (non-manifold boundary); leave it to the check
      nleft = nleft + 1;
      continue
    end
    if numel(loop) > maxloop
      nleft = nleft + 1;
      continue
    end
    % fan around the centroid; the new triangle reverses the boundary edge
    % (a->b becomes b->a) so that its orientation matches the existing face
    pos(end + 1, :) = mean(pos(loop, :), 1); %#ok<AGROW>
    c = size(pos, 1);
    for j = 1:numel(loop)
      a = loop(j);
      b = loop(mod(j, numel(loop)) + 1);
      tri(end + 1, :) = [b a c]; %#ok<AGROW>
    end
    nclosed = nclosed + 1;
  end

  bnd(ii).pos = pos;
  bnd(ii).tri = tri;
  if nclosed > 0 || nleft > 0
    fprintf(['close_mesh_holes: mesh %d had %d boundary edges, closed %d hole(s)', ...
             ', %d left open\n'], ii, size(boundary, 1), nclosed, nleft);
  end
end
end %close_mesh_holes
