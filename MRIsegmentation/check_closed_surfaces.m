function check_closed_surfaces(bnd, label, expected_ml)
% CHECK_CLOSED_SURFACES  Fail when a surface mesh is not watertight.
%
% A boundary surface that is not closed cannot be used for BEM, for deciding
% whether a point lies inside the head, or for placing optodes on the scalp,
% and the defect is easy to miss: cgalsurf returns such meshes without any
% warning. Called after close_mesh_holes, so anything reported here is a real
% problem rather than a small meshing glitch.
%
% Checked per mesh: no boundary edges, no non-manifold edges, outward
% orientation (positive enclosed volume) and, when expected_ml is given, that
% the mesh encloses the whole compartment rather than just a layer of it - a
% double-walled shell is closed and consistently wound and would otherwise
% pass unnoticed.
%
% Input:   bnd         <struct array> meshes with .pos and .tri
%          label       <string> name used in the message (optional)
%          expected_ml <vector> compartment volume per mesh in ml (optional)
%
% (c) Eike Middell, October 2026

if nargin < 2
  label = 'surface';
end

if nargin < 3
  expected_ml = [];
end

bad = {};
for ii = 1:numel(bnd)
  tri = bnd(ii).tri;
  pos = bnd(ii).pos;
  enclosed = sum(dot(pos(tri(:, 1), :), cross(pos(tri(:, 2), :), pos(tri(:, 3), :), 2), 2)) / 6000;
  edges = sort([tri(:, [1 2]); tri(:, [2 3]); tri(:, [3 1])], 2);
  [~, ~, idx] = unique(edges, 'rows');
  counts = accumarray(idx, 1);
  nopen = nnz(counts == 1);
  nnonmanifold = nnz(counts > 2);
  fprintf(['check_closed_surfaces: %s %d: %d vertices, %d faces, %d boundary edges, ' ...
           '%d non-manifold edges, encloses %.0f ml\n'], ...
          label, ii, size(pos, 1), size(tri, 1), nopen, nnonmanifold, enclosed);
  if nopen > 0 || nnonmanifold > 0
    bad{end + 1} = sprintf('%d: %d boundary, %d non-manifold edges', ii, nopen, nnonmanifold); %#ok<AGROW>
  elseif enclosed <= 0
    bad{end + 1} = sprintf('%d: inside out (encloses %.0f ml)', ii, enclosed); %#ok<AGROW>
  elseif ~isempty(expected_ml) && enclosed < 0.8 * expected_ml(ii)
    bad{end + 1} = sprintf(['%d: encloses only %.0f of %.0f ml, i.e. a shell around the ' ...
                            'compartment rather than the compartment itself'], ...
                           ii, enclosed, expected_ml(ii)); %#ok<AGROW>
  end
end

if ~isempty(bad)
  error('check_closed_surfaces:open', ...
        ['unusable %s mesh - %s. Small holes are closed by close_mesh_holes and the ' ...
         'orientation by orient_outward, so what remains points at the meshing itself: ' ...
         'too few vertices (raise num_vertices), or a compartment that is not solid.'], ...
        label, strjoin(bad, '; '));
end
end %check_closed_surfaces
