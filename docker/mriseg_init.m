function mriseg_init()
% MRISEG_INIT  Put all pipeline dependencies on the MATLAB path.
%
% Called by the `mriseg` wrapper before every command inside the container.
% Changes the working directory to $MRISEG_HOME because the pipeline scripts
% resolve their resources (Huang_et_al_2013/eTPM.nii, fieldtrip/external/spm12/
% toolbox/cat12, Schaefer2018_Parcellations, HArtMuT) relative to pwd.
%
% Deliberately does NOT genpath the whole FieldTrip tree: its compat/ and
% external/ folders shadow MATLAB builtins and break SPM12 (see README).

home = getenv('MRISEG_HOME');
if isempty(home)
  home = '/opt/MRIsegmentation';
end
cd(home);

addpath(home);

% FieldTrip: ft_defaults adds the required subfolders itself
addpath(fullfile(home, 'fieldtrip'));
ft_defaults;
% helpers the pipeline calls directly from FieldTrip's external/ folder
addpath(fullfile(home, 'fieldtrip', 'external', 'openmeeg'));    % om_save_tri
addpath(fullfile(home, 'fieldtrip', 'external', 'freesurfer'));  % read_curv

% SPM12 + CAT12 (full SPM12 replaces FieldTrip's trimmed copy in the image)
spm12 = fullfile(home, 'fieldtrip', 'external', 'spm12');
addpath(spm12);
addpath(fullfile(spm12, 'toolbox', 'cat12'));
spm('defaults', 'fmri');
spm_get_defaults('cmdline', true);

addpath(fullfile(home, 'iso2mesh'));
addpath(fullfile(home, 'Huang_et_al_2013'));
addpath(fullfile(home, 'NIfTI_tools'));
addpath(fullfile(home, 'weighted-SPHARM'));

% last, so the pipeline's own folder shadows everything above
addpath(fullfile(home, 'MRIsegmentation'));

fprintf('MRIsegmentation environment ready (%s, MATLAB %s)\n', home, version);
end
