function start_cat(NiiFile, atlas, Template)
%
% Input:        NiiFile <string> fullpath to T1 MRI image to be segmented.
%               atlas <string> fullpath to cortical atlas used by CAT12.
%               Template <string> fullpath to prior probability distribution
%                                 template used for segmentation. 
%
% This function calls the CAT12 segmentation routine. It does the following:
% 1. Sets CAT12 segmentation settings.
% 2. Defines batch jobs for parallel computing. 
% 3. Start batch jobs.
%
% Example:
% read_cat('/usr/local/cat12', '/tmp/head/1/T1.nii', 2000);
% start_cat('/tmp/head/1/T1.nii', '/usr/local/CAT12/atlases_surfaces/...
%                                  lh.aparc_DK40.freesurfer.annot', 
%                                 'eTPM.nii');
% 
% (c) Nils Harmening, May 2020
% Neurotechnology group, Technische Universität Berlin, Germany

if nargin <1 || isempty(NiiFile) %no files
    NiiFile = spm_select(inf,'image','Select images for new segment');
end

%% Create SPM batch
matlabbatch{1}.spm.tools.cat.estwrite.data = {[NiiFile ',1']};
matlabbatch{1}.spm.tools.cat.estwrite.nproc = 0;
matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm = {atlas};
matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg = 'mni';
matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.opts.accstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.APP = 1070;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.gcutstr = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.dartel.darteltpm = {Template};
matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox = 1.5;

%% Version 12.7
matlabbatch{1}.spm.tools.cat.estwrite.extopts.restypes.fixed = [1 0.1];
matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.noROI = struct([]);
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40 = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod = 0; %1
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod = 0; %1
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.labelnative = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped = 0; %1
matlabbatch{1}.spm.tools.cat.estwrite.output.jacobianwarped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.warps = [0 0];

% Run SPM batch
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);

end %start_cat
