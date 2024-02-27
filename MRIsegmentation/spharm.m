function spharm(input_img, L)
%
% Input:        input_img <string> fullpath to T1 MRI image to be segmented
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
% spharm('/tmp/head/1/T1.nii')
% 
% (c) Nils Harmening, February 2024
% Neurotechnology group, Technische Universität Berlin, Germany

[dirname, base_filename, ext] = fileparts(input_img);


%% Load hemisphere surfaces and their corresponding sphere meshes
load(fullfile(dirname, 'ctf', 'lh_surf.mat'))
bnd_lh.vertices = bnd_lh.pos;
bnd_lh.faces = bnd_lh.tri;
load(fullfile(dirname, 'ctf', 'lh_sphere.mat'))
lh_sphere.vertices = lh_sphere.pos;
lh_sphere.faces = lh_sphere.tri;

load(fullfile(dirname, 'ctf', 'rh_surf.mat'))
bnd_rh.vertices = bnd_rh.pos;
bnd_rh.faces = bnd_rh.tri;
load(fullfile(dirname, 'ctf', 'rh_sphere.mat'))
rh_sphere.vertices = rh_sphere.pos;
rh_sphere.faces = rh_sphere.tri;



%% Get the weighted-SPHARM representation. In order to smooth with up to 78th degree weighted-SPHARM with bandwidth 0.0001 (Figure 1), run 
sigma=0.0001;
[lh_surf_smooth, lh_fourier_coeff]=SPHARMsmooth2(bnd_lh,lh_sphere,L,sigma);
[rh_surf_smooth, rh_fourier_coeff]=SPHARMsmooth2(bnd_rh,rh_sphere,L,sigma);


%% Extracting Spherical harmonic coefficients
k = L; %78;
lh_fourier_vector=SPHARMvectorize(lh_fourier_coeff,k);
save(fullfile(dirname, 'ctf', 'lh_fourier_coeff.mat'), 'lh_fourier_coeff', '-v7');
save(fullfile(dirname, 'ctf', 'lh_fourier_vector.mat'), 'lh_fourier_vector', '-v7');

rh_fourier_vector=SPHARMvectorize(rh_fourier_coeff,k);
save(fullfile(dirname, 'ctf', 'rh_fourier_coeff.mat'), 'rh_fourier_coeff', '-v7');
save(fullfile(dirname, 'ctf', 'rh_fourier_vector.mat'), 'rh_fourier_vector', '-v7');

end %spharm
