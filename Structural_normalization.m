% This code coregister and normalizes the MVPA results from naive space to
% MNI standard space. It is based on the steps in, except the very last
% smoothing step:
% https://www.youtube.com/watch?v=zSqBoB1GrDk

% 1) Coregisteration: Ref: meanaV0001, souce: V0001 anatomical image, previosuly oreinted AC, Output: WV0001
% 2) Segmentation: Save biase corrected, Forward & Inverse
% 3) Normalization (write): Def: y_rV0001, images to be realigned: MVPA outputs

clear
spm('defaults','fmri');
spm_jobman('initcfg');

AnatFolder='C:\Users\user\Documents\Reyhaneh\EEG+fMRI\EO\MRI_data\Anat';
FuncFolder='C:\Users\user\Documents\Reyhaneh\EEG+fMRI\EO\MRI_data\Func';
MVPAFolder='C:\Users\user\Documents\Reyhaneh\EEG+fMRI\EO\MRI_data\Analysis\MVPA\SpotlightR0';
Subjects = {'AB_04152014' 'AM_05082014' 'JK_04252014' 'KK_04152014' 'MB_04252014' 'MG_04152014' 'ML_05122014' 'MM_05132014' 'NB_05012014' ...
    'TD_05012014' 'CJ_06052014'  'CS_04162014' 'DH_05142014' 'NJM_05142014' 'PD_04162014' ...
'RS_04162014' 'SA_05162014' 'TP_05162014' 'DS_12132013' ...
'ET_12102013' 'RS_12102013' 'TH_12172013'};



%% Coregisteration
COREG=0;
if (COREG==1)
Coregisteration_batchTemplate=load('C:\Users\user\Documents\Reyhaneh\EEG+fMRI\Codes\Matlab\Coregisteration_batch_template.mat');
numPar=1;
for SubIDInd=1:length(Subjects)
    SubID=Subjects{SubIDInd};
    disp(SubID)
    RefImg = [FuncFolder,'\',SubID,'\run01\meanaV0001.img,1']; %representation image, mean funct
    SourceImg = [AnatFolder,'\',SubID,'\V0001.img,1']; %Anatomical, AC reoriented
    
    matlabbatch{numPar}=Coregisteration_batchTemplate.matlabbatch{1};
    matlabbatch{numPar}.spm.spatial.coreg.estwrite.ref{1} = RefImg;
    matlabbatch{numPar}.spm.spatial.coreg.estwrite.source{1} = SourceImg;
    
%    numPar=numPar+1;
    spm_jobman('run',matlabbatch);
end

end

%% Segmentation 
SEG=1;
if (SEG==1)
Segmentation_batchTemplate=load('C:\Users\user\Documents\Reyhaneh\EEG+fMRI\Codes\Matlab\Segmentation_batch_template.mat');

numPar=1;
for SubIDInd=1 :length(Subjects)
    SubID=Subjects{SubIDInd};
    disp(SubID)
   
    RealignedStructuralImage=[AnatFolder,'\',SubID,'\rV0001.img,1']; %Realigned Anatomical, AC reoriented
    
    matlabbatch{numPar}=Segmentation_batchTemplate.matlabbatch{1};
    matlabbatch{numPar}.spm.spatial.preproc.channel.vols{1}=RealignedStructuralImage;

%    numPar=numPar+1;
    spm_jobman('run',matlabbatch);
end
%spm_jobman('run',matlabbatch);
end

%% Normalization
NORM=1;
if (NORM==1)
Normalization_batchTemplate=load('C:\Users\user\Documents\Reyhaneh\EEG+fMRI\Codes\Matlab\Normalization_batch_template.mat');

numPar=1;
for SubIDInd=1 :length(Subjects)
    SubID=Subjects{SubIDInd};
    disp(SubID)
   
    DefImage=[AnatFolder,'\',SubID,'\y_rV0001.nii']; % 
    mage_2_warp1=[FuncFolder,'\',SubID,'\run01\meanaV0001.img,1']; % image in naive space to be wrapped into MNI space (meanaV0001 to check Normalization accuracy)
    %Unzip MVPA output image
    Image_2_warp2gz=[MVPAFolder,'\MVPA_',SubID,'.nii.gz'];
    gunzip(Image_2_warp2gz);
    Image_2_warp2=[MVPAFolder,'\MVPA_',SubID,'.nii']; % image in naive space to be wrapped into MNI space (MVPA image)
    
    matlabbatch{numPar}=Normalization_batchTemplate.matlabbatch{1};
    matlabbatch{numPar}.spm.spatial.normalise.write.subj.def{1}=DefImage;
    matlabbatch{numPar}.spm.spatial.normalise.write.subj.resample{1}=mage_2_warp1;
    matlabbatch{numPar}.spm.spatial.normalise.write.subj.resample{2}=Image_2_warp2;
    
%    numPar=numPar+1;
    spm_jobman('run',matlabbatch);
    
end
%spm_jobman('run',matlabbatch);
end