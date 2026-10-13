function bnd = orient_outward(bnd)
% ORIENT_OUTWARD  Make the triangles of closed surfaces face outwards.
%
% iso2mesh/cgalsurf does not guarantee a particular orientation, and FieldTrip
% does not fix it: the skull surface of a run can come out inside-out while
% the others are fine. Normals pointing the wrong way break BEM solvers and
% any inside/outside test. For a closed, consistently wound mesh the sign of
% the enclosed volume tells the orientation, so flipping all triangles of the
% meshes with negative volume is enough.
%
% (c) Eike Middell, October 2026

for ii = 1:numel(bnd)
  pos = bnd(ii).pos;
  tri = bnd(ii).tri;
  v = sum(dot(pos(tri(:, 1), :), cross(pos(tri(:, 2), :), pos(tri(:, 3), :), 2), 2)) / 6;
  if v < 0
    bnd(ii).tri = tri(:, [1 3 2]);
    fprintf('orient_outward: mesh %d was inside out (volume %.0f ml), flipped\n', ii, v / 1000);
  end
end
end %orient_outward
