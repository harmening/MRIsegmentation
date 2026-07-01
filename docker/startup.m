% MRIsegmentation MATLAB startup script
disp('Loading MRIsegmentation environment...');

% Change to MRIsegmentation directory (scripts expect this as working dir)
cd('/opt/MRIsegmentation');

% Add paths
addpath('/opt/MRIsegmentation');
addpath(genpath('/opt/MRIsegmentation/MRIsegmentation'));
addpath(genpath('/opt/MRIsegmentation/fieldtrip'));
addpath(genpath('/opt/MRIsegmentation/iso2mesh'));
addpath(genpath('/opt/MRIsegmentation/Huang_et_al_2013'));
addpath(genpath('/opt/MRIsegmentation/NIfTI_tools'));
addpath(genpath('/opt/MRIsegmentation/weighted-SPHARM'));
addpath(genpath('/opt/MRIsegmentation/Schaefer2018_Parcellations'));

% Initialize FieldTrip (adds necessary subfolders)
ft_defaults;

% Configure iso2mesh
setenv('ISO2MESH_BIN', '/opt/MRIsegmentation/iso2mesh/bin');

disp('MRIsegmentation environment loaded.');
