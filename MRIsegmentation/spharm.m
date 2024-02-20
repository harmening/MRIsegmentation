function spharm(input_img, spharm_dir)
%
% Input:        input_img <string> fullpath to T1 MRI image to be segmented
% Input:        spharm_dir <string> fullpath to spherical harmonics (will be
%                                   created if not existing)
% 
% This function constructs the weighted-SPHARM representation from the cortex
% surface (output of read_cat.m). It does the follwoing:
% 1. Reads the cortex surface.
% 2. Computes spherical harmonics Y_lm from the mesh of a unit sphere
%    (unitsphere.mat). The weighted-SPHARM representation is built upon
%    this mesh.
% 3. Get the weighted-SPHARM representation for input cortex surface. 
% 4. Extract spherical harmonic coefficients and stores them as vector
%    representation.
%
% Example:
% spharm('/tmp/head/1/T1.nii', '~/MRIsegmentation/weighted-SPHARM/sph/sph');
% 
% (c) Nils Harmening, February 2024
% Neurotechnology group, Technische Universität Berlin, Germany

[dirname, base_filename, ext] = fileparts(input_img);

%% Construct SPHARM
if ~numel(strcat(spharm_dir, '0.mat'))
  SPHARMconstruct(spharm_dir, 85);
end


%% Load cortex surface
load(fullfile(dirname, 'ctf', 'cortex_surf40962.mat'));
coord = cortex_surf.pos';
tri = cortex_surf.tri;


%% Get the weighted-SPHARM representation. In order to smooth with up to 78th degree weighted-SPHARM with bandwidth 0.0001 (Figure 1), run 
L=78;
sigma=0.0001;
[coord_smooth,fourier_coeff] = SPHARMsmooth(spharm_dir, coord,L,sigma);


%% Extracting Spherical harmonic coefficients
k = L; %50;
fourier_vector=SPHARMvectorize(fourier_coeff,k)

save(fullfile(dirname, 'ctf', 'fourier_coeff.mat'), 'fourier_coeff', '-v7');
save(fullfile(dirname, 'ctf', 'fourier_vector.mat'), 'fourier_vector', '-v7');

end %spharm
