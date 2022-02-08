function [new_pos] = applyDeformation(deffilename, pos)
%
% Input:        deffilename <string> fullpath to the SPM deformation file
%               pos <nx3 float> point coordinates to be transformed.
% Output:       new_pos <nx3 float> transformed point coordinates.
%
% This function applies the nonlinear SPM voxelbased transformation as stored
% in deffilename to points specified in pos. It does the following:
% 1. Loads the nonlinear transformation matrix and inverts it.
% 2. Interpolates the mapping between the values.
% 3. Transformes the input points. 
%
% Example:
% applyDeformation('y_example_T1.nii', points); 
% 
% (c) Nils Harmening, Daniel Miklody, May 2020
% Neurotechnology group, Technische Universität Berlin, Germany

%% Load deformation as produced by spm8 new_segment
Vol=spm_vol([repmat(deffilename,3,1),[',1,1';',1,2';',1,3']]);

%% Inverse transform of template
M_initial.transform = Vol(1).mat;
Mtemplate=inv(M_initial.transform);
pos=ft_warp_apply(Mtemplate, pos, 'homogeneous');

Def=zeros([Vol(1).dim 3]);
Def(:,:,:,1)= spm_load_float(Vol(1));
Def(:,:,:,2)= spm_load_float(Vol(2));
Def(:,:,:,3)= spm_load_float(Vol(3));

%% Create grid
[X, Y, Z]=ndgrid(1:Vol(1).dim(1),1:Vol(1).dim(2),1:Vol(1).dim(3));

%% Construct interpolants
V = Def(:,:,:,1);
%XInterp = griddedInterpolant(X,Y,Z,V,'nearest');
XInterp = griddedInterpolant(X,Y,Z,V);
V = Def(:,:,:,2);
%YInterp = griddedInterpolant(X,Y,Z,V,'nearest');
YInterp = griddedInterpolant(X,Y,Z,V);
V = Def(:,:,:,3);
%ZInterp = griddedInterpolant(X,Y,Z,V,'nearest');
ZInterp = griddedInterpolant(X,Y,Z,V);

%% Apply interpolants to pos
new_pos=[XInterp(pos) YInterp(pos) ZInterp(pos)];
% new_pos=ft_warp_apply(M,new_pos,'homogeneous');

end %applyDeformation
