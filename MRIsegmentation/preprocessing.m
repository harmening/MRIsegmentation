function preprocessing(mri_lst)
%
% Input:        mri_lst <list of strings> list of fullpaths to MRI images to be
%                                         preprocessed
%
% This function converts the MRI scan to nifti-fileformat, orients the input
% MRIs to the ACPC corrdinate system and creates an average MRI image (useful
% for multiple scans of the same subject). It does the following:
% 1. Converts the input images to nifti, if necessary, using freesurfers
%    mri_convert.
% 2. Calls ATRA, which detects anatomical landmarks like AC, PC, VSPS and
%    reorientates the MRI into the ACPC coordinate system (origin in AC,
%    orientation RAS, individual geometry i.e. not normalized). 
%(3.)For a correction of missing values in the MRI scans, see python-based
%    imputation parallel.py
%
% Example:
% preprocessing(['/tmp/head/1/T1_1.img', '/tmp/head/1/T1_2.img']); 
% Results are stored in /tmp/head/1/ and include the nifti-images T1_1.nii,
% T1_2.nii, their ACPC-aligned copy T1_1_RAS.nii, T1_2_RAS.nii and their
% average atra_avg_RAS.nii.
% 
% (c) Nils Harmening, May 2020
% Neurotechnology group, Technische Universität Berlin, Germany

% Conversion to nifti with mri_convert of freesurfer 
%for i=numel(mri_lst)
%    mri_img = mri_lst{i};
%    [filepath, filename, ext] = fileparts(mri_img);
%    if strcmp(ext, '.nii')
%        mri_nii = fullfile(filepath, strcat(filename, '.nii'));
%        system(['mri_convert ', mri_img, ' ', mri_nii]);
%        mri_lst(i) = mri_nii;
%    end
%end

% Correction for missing values (NOT TESTED!!!!!) 
%for i=numel(mri_lst)
%  mri_img = mri_lst{i};
%  mri = ft_read_mri(mri_img);
%  mri(isnan(mri))=0;
%  save_nii(mri_img, mri);
%end

% Prepare atra orientation detection + correction
[filepath, filename, ext] = fileparts(mri_lst{1});
imagelist = fullfile(filepath, 'imagelist');
fileID = fopen(imagelist, 'w');
for i=numel(mri_lst)
    fprintf(fileID, strcat(mri_lst{i},'\n'));
end
fclose(fileID);

% Run atra
options = '-v -orient RAS';
system(['cd ', filepath, ';', 'atra ', options, ' -i ', imagelist]);

% atra clean up
output = {strcat(filename, '*'), imagelist, 'atra.txt', 'atra_avg_RAS.nii'};
for i=numel(mri_lst)
    output{end+1} = mri_lst(i);
end 
for file = output
    try
      system(['mv ', file, ' ', filepath]);
    catch
      disp('Errors during ATRA clean up');
    end
end

end %preprocessing
