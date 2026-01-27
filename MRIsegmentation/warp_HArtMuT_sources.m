function warp_HArtMuT_sources_iso2mesh( ...
    srcScalpFile, srcSkullFile, ...
    tgtScalpFile, tgtSkullFile, ...
    hartmutMatIn, hartmutMatOut, ...
    mean_pnt)
%WARP_HARTMUT_SOURCES_ISO2MESH
% Warp HArtMuT artefact source positions from a source head (e.g., NYHead)
% to an individual target head using radial depth transfer between skull
% and scalp surfaces along rays from a fixed mean point.
%
% - Supports mesh input formats: .stl and OpenMEEG .tri (written by om_save_tri.m)
% - Uses iso2mesh raytrace() for ray/triangle intersections.
%
% Inputs:
%   srcScalpFile  - source scalp mesh (.stl or .tri)
%   srcSkullFile  - source skull mesh (.stl or .tri)
%   tgtScalpFile  - target scalp mesh (.stl or .tri) (neck-extended scalp)
%   tgtSkullFile  - target skull mesh (.stl or .tri)
%   hartmutMatIn  - input .mat containing HArtMuT.artefactmodel.pos, etc.
%   hartmutMatOut - output .mat with individualized positions
%   mean_pnt      - 1x3 ray origin (default [0 -10 0])
%
% Requirements:
%   - iso2mesh on MATLAB path (raytrace.m must be resolvable)
%   - stlread on path (MATLAB built-in on newer versions; otherwise add one)
%
% Example:
%   mean_pnt = [0 -10 0];
%   warp_HArtMuT_sources_iso2mesh( ...
%     './NYhead/scalp.tri', './NYhead/skull.tri', ...
%     './Colin27/extended_scalp.tri', './Colin27/skull.tri', ...
%     '../HArtMuTmodels/HArtMuT_NYhead_small.mat', ...
%     './Colin27/HArtMuT_Colin27.mat', ...
%     mean_pnt);

if nargin < 7 || isempty(mean_pnt)
    mean_pnt = [0 -10 0];
end
mean_pnt = double(mean_pnt(:).'); % 1x3

if exist('raytrace', 'file') ~= 2
    error('iso2mesh raytrace.m not found on your MATLAB path.');
end

% --- Load meshes
[srcScalpV, srcScalpF] = read_surface_mesh(srcScalpFile);
[srcSkullV, srcSkullF] = read_surface_mesh(srcSkullFile);
[tgtScalpV, tgtScalpF] = read_surface_mesh(tgtScalpFile);
[tgtSkullV, tgtSkullF] = read_surface_mesh(tgtSkullFile);

% --- Load HArtMuT model
S = load(hartmutMatIn);
if ~isfield(S, 'HArtMuT') || ~isfield(S.HArtMuT, 'artefactmodel')
    error('Input MAT must contain HArtMuT.artefactmodel');
end
artefactmodel = S.HArtMuT.artefactmodel;

if ~isfield(artefactmodel, 'pos')
    error('HArtMuT.artefactmodel.pos not found.');
end
sources = double(artefactmodel.pos);
if size(sources,2) ~= 3
    error('artefactmodel.pos must be Nx3');
end

% --- Warp points
new_sources = zeros(size(sources));
n = size(sources,1);

for i = 1:n
    p = sources(i,:);
    d = p - mean_pnt;
    dn = norm(d);
    if dn < 1e-9
        new_sources(i,:) = p;
        continue;
    end
    dir = d / dn;
    t_p = dn;

    % Intersections in source head
    tSk_src = first_positive_intersection_iso2mesh(mean_pnt, dir, srcSkullV, srcSkullF);
    tSc_src = first_positive_intersection_iso2mesh(mean_pnt, dir, srcScalpV, srcScalpF);

    % Intersections in target head
    tSk_tgt = first_positive_intersection_iso2mesh(mean_pnt, dir, tgtSkullV, tgtSkullF);
    tSc_tgt = first_positive_intersection_iso2mesh(mean_pnt, dir, tgtScalpV, tgtScalpF);

    if any(isnan([tSk_src, tSc_src, tSk_tgt, tSc_tgt]))
        % Fallback: keep original
        new_sources(i,:) = p;
        continue;
    end

    % Ensure skull is inside scalp along the ray (swap if needed)
    if tSk_src > tSc_src, [tSk_src, tSc_src] = deal(tSc_src, tSk_src); end
    if tSk_tgt > tSc_tgt, [tSk_tgt, tSc_tgt] = deal(tSc_tgt, tSk_tgt); end

    denom = (tSc_src - tSk_src);
    if abs(denom) < 1e-9
        new_sources(i,:) = p;
        continue;
    end

    % Relative depth between skull and scalp in source
    alpha = (t_p - tSk_src) / denom;

    % Apply same relative depth in target
    t_new = tSk_tgt + alpha * (tSc_tgt - tSk_tgt);
    new_sources(i,:) = mean_pnt + t_new * dir;
end

% --- Save individualized model
individual_artefactmodel = artefactmodel;
individual_artefactmodel.pos = new_sources;
fields = {'leadfield', 'solver', 'modeltype'};
individual_artefactmodel = rmfield(individual_artefactmodel, fields);

artefactmodel = individual_artefactmodel;

save(hartmutMatOut, 'artefactmodel');
fprintf('Wrote individualized HArtMuT model: %s\n', hartmutMatOut);

end

% ======================================================================
% Helper: closest positive ray/mesh intersection using iso2mesh raytrace
% ======================================================================
function tMin = first_positive_intersection_iso2mesh(orig, dir, V, F)
[t, u, v, idx] = raytrace(orig, dir, V, F); %#ok<ASGLU>
if isempty(idx)
    tMin = NaN;
    return;
end
tt = t(idx);
tt = tt(isfinite(tt) & tt > 1e-6);
if isempty(tt)
    tMin = NaN;
else
    tMin = min(tt);
end
end

% ======================================================================
% Mesh readers: .stl and OpenMEEG .tri
% ======================================================================
function [V, F] = read_surface_mesh(filename)
[~,~,ext] = fileparts(filename);
ext = lower(ext);

switch ext
    case '.stl'
        if exist('stlread','file') ~= 2
            error('stlread not found on MATLAB path.');
        end
        out = stlread(filename);
        if isa(out, 'triangulation')
            F = out.ConnectivityList;
            V = out.Points;
        else
            % Some stlread variants return [F,V]
            try
                [F,V] = stlread(filename);
            catch
                error('stlread output format not recognized.');
            end
        end

    case '.tri'
        [V, F] = read_tri(filename);

    otherwise
        error('Unsupported mesh format: %s (use .stl or .tri)', ext);
end

V = double(V);
F = double(F);
end

function [V, F] = read_tri(filename)
% Inverse of om_save_tri.m
fid = fopen(filename, 'r');
if fid < 0
    error('Cannot open TRI file: %s', filename);
end
cleanup = onCleanup(@() fclose(fid));

% number of points
line = fgetl(fid);
tmp  = sscanf(line, '- %d');
npoints = tmp(1);

% points + normals (6 columns per vertex)
data = fscanf(fid, '%f', [6 npoints]);
data = data.';
V = data(:,1:3); % normals ignored

% number of faces
line = fgetl(fid);
if length(line) < 3 %empty
    line = fgetl(fid);
end
tmp  = sscanf(line, '- %d %d %d');
nfaces = tmp(1);

% faces (0-based in file)
faces = fscanf(fid, '%d', [3 nfaces]);
faces = faces.' + 1;
F = faces;
end
