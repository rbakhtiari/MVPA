% SPM12 batch script for performing 1st level analysis
% for MVPA analysis. beta will be calculated for 'a' files (naive space), and
% each stimuli will be associated with a seperate beta, vs all stimuli type
% be associated with a seperate beta as in traditional univariate apporach
% for each session a seperate SPM is generated as the ra files are in naive
% space of that subject session.

%% Initial Setting
clear
spm('defaults','fmri');
spm_jobman('initcfg');

%--- BASIC SETUP ---%

% subjects = {'AM_05082014'}; %run order 1-5
%subjects = {'AB_04152014' 'JK_04252014' 'KK_04152014' 'MB_04252014' 'MG_04152014' 'ML_05122014' 'MM_05132014' 'NB_05012014' 'TD_05012014' 'CJ_06052014' 'CS_04162014' 'DH_05142014' 'NJM_05142014' 'PD_04162014' 'RS_04162014' 'SA_05162014' 'TP_05162014'}; %5 runs
%subjects = {'AB_04152014' 'CJ_06052014' 'CS_04162014' 'DH_05142014' 'JK_04252014' 'MM_05132014' 'SA_05162014' 'TD_05012014' 'TP_05162014'};
%subjects = {'DS_12132013' 'ET_12102013' 'RS_12102013' 'TH_12172013'}; %5 runs
%subjects = {'TP_05162014'}; %5 runs
%subjects = {'DS_12132013' 'ET_12102013' 'RS_12102013' 'TH_12172013'}; %8 runs

subjects = {'AB_04152014' 'AM_05082014' 'JK_04252014' 'KK_04152014' 'MB_04252014' 'MG_04152014' 'ML_05122014' 'MM_05132014' 'NB_05012014' ...
    'TD_05012014' 'CJ_06052014'  'CS_04162014' 'DH_05142014' 'NJM_05142014' 'PD_04162014' ...
'RS_04162014' 'SA_05162014' 'TP_05162014' 'DS_12132013' ...
'ET_12102013' 'RS_12102013' 'TH_12172013'};

analysis = '3cond_dist2';
CondTypeName={'Neu','Neg','Pos','Targ'};

FuncFolder='C:\Users\user\Documents\Reyhaneh\EEG+fMRI\EO\MRI_data\Func\';
rootdir = 'C:\Users\user\Documents\Reyhaneh\EEG+fMRI\EO\MRI_data\Analysis\SPM';
FirstLevelFolder='C:\Users\user\Documents\Reyhaneh\EEG+fMRI\EO\MRI_data\FirstLevel\';

 
%% --- Realignment---%
Realignment=0;
if (Realignment==1)
    clear matlabbatch
    % All sessions are entered, registered to the first image, for reslicing all images+mean image is resliced
    reslilce_batchTemplate=load('C:\Users\user\Documents\Reyhaneh\EEG+fMRI\Codes\Matlab\Realign_EstReslice_AllSessions.mat');
   % reslilce_batchTemplate=load('C:\Users\user\Documents\Reyhaneh\EEG+fMRI\Codes\Matlab\Realign_scr-old.mat');
    numPar=1;
    for iSubjID = 1:length(subjects)
        SubjID=subjects{iSubjID};
        
        SPM=load([FirstLevelFolder,SubjID,'\3cond_distr2\SPM.mat']);
        matlabbatch{numPar}=reslilce_batchTemplate.matlabbatch{1};
    
        nRuns=length(SPM.SPM.Sess);
        disp([SubjID, ': ', num2str(nRuns)])
        for iSess =1:nRuns
            disp(['Session: ',num2str(iSess)])
    
            % Functional Image:
            FirstSessImgRowID=SPM.SPM.Sess(iSess).row(1);
            LastSessImgRowID=SPM.SPM.Sess(iSess).row(end);
    
            FirstImg=SPM.SPM.xY.P(FirstSessImgRowID,:);
            RunID=['run',num2str(iSess,'%02.f')];
            BaseEnd=strfind(FirstImg,'Func')+4;
            ImageStart=strfind(FirstImg,'swaV')-1;
            SubjRunFolder=FirstImg(BaseEnd+1:ImageStart);
            FuncName={};
            cntfName=1;
            for iImage=FirstSessImgRowID:LastSessImgRowID
                
                ImgID=SPM.SPM.xY.P(iImage,ImageStart+3:end);
                MyImage=[FuncFolder,SubjRunFolder,ImgID];
                FuncName{cntfName}=MyImage;
                if ((iImage==FirstSessImgRowID) || (iImage==LastSessImgRowID))
                    if ~exist(MyImage)
                        disp('********************************************')
                        disp('  Error! Func file does not exist!!')
                    end
                end
                 
                cntfName=cntfName+1;
            end
            FuncName=FuncName';
            disp(['   Total functional file: ',num2str(cntfName-1)])
            matlabbatch{numPar}.spm.spatial.realign.estwrite.data{iSess} = FuncName;
            
        end
     %   numPar=numPar+1;
        spm_jobman('run',matlabbatch);
    end
end
%spm_jobman('run',matlabbatch);

%% --- ART---%
doART=0;

ARTrootdir = 'C:\Users\user\Documents\Reyhaneh\EEG+fMRI\EO\MRI_data\Analysis\ART';
                                                                                                                                                                               
%%%%%%%%%%%% ART PARAMETERS (edit to desired values) %%%%%%%%%%%%
global_mean=1;                % global mean type (1: Standard 2: User-defined Mask)
motion_file_type=0;           % motion file type (0: SPM .txt file 1: FSL .par file 2:Siemens .txt file)
global_threshold=9.0;         % threshold for outlier detection based on global signal
motion_threshold=2.0;         % threshold for outlier detection based on motion estimates
% use_diff_motion=1;            % 1: uses scan-to-scan motion to determine outliers; 0: uses absolute motion
% use_diff_global=1;            % 1: uses scan-to-scan global signal change to determine outliers; 0: uses absolute global signal values
% use_norms=1;                  % 1: uses composite motion measure (largest voxel movement) to determine outliers; 0: uses raw motion measures (translation/rotation parameters) 
% mask_file=[];                 % set to user-defined mask file(s) for global signal estimation (if global_mean is set to 2) 

 
if (doART==1)
    for iSubjID = 1:length(subjects)
        SubjID=subjects{iSubjID};
        
        SPM=load([FirstLevelFolder,SubjID,'\3cond_distr2\SPM.mat']);
    
        nRuns=length(SPM.SPM.Sess);
        disp([SubjID, ': ', num2str(nRuns)])

        cfgfile=fullfile(ARTrootdir,['art_config_',SubjID,'.cfg']);
        fid=fopen(cfgfile,'wt');
        %[filepath,filename,fileext]=fileparts(deblank(files(n1,:)));
        
        fprintf(fid,'# Users can edit this file and use\n');
        fprintf(fid,'#   art(''sess_file'',''%s'');\n',cfgfile);
        fprintf(fid,'# to launch art using this configuration\n');
        
        fprintf(fid,'sessions: %d\n',nRuns);
        fprintf(fid,'global_mean: %d\n',global_mean);
        fprintf(fid,'global_threshold: %f\n',global_threshold);
        fprintf(fid,'motion_threshold: %f\n',motion_threshold);
        fprintf(fid,'motion_file_type: 0\n');
        fprintf(fid,'motion_fname_from_image_fname: 0\n');
%         fprintf(fid,'use_diff_motion: %d\n',use_diff_motion);
%         fprintf(fid,'use_diff_global: %d\n',use_diff_global);
%         fprintf(fid,'use_norms: %d\n',use_norms);
%         fprintf(fid,'output_dir: %s\n',art_fileparts(files(n1,:)));
%         if ~isempty(mask_file),fprintf(fid,'mask_file: %s\n',deblank(mask_file(n1,:)));end
        image_dir=[FuncFolder,SubjID];
        fprintf(fid,'image_dir: %s\n', image_dir); % functional and movement data folders (comment these lines if functional/movement filenames below contain full path information)
        fprintf(fid,'motion_dir: %s\n',image_dir);

        fprintf(fid,'end\n');
        fprintf(fid,'\n');
        
        
        for iSess=1:nRuns
            fprintf(fid,'session %d image run0%d/raV0???.img\n',iSess,iSess); % 'session 1 image run01/av0???.img' 
        end 
        fprintf(fid,'\n');

        for iSess=1:nRuns
            fprintf(fid,'session %d motion run0%d/rp_aV0001.txt\n',iSess,iSess); % 'session 1 motion run01/rp_aV0001.txt '
        end
        fprintf(fid,'\n');

        fprintf(fid,'end\n');
        fclose(fid);
    
        art('sess_file',cfgfile);
    end

    
end

%% --- MODEL SPECIFICATION, and contrast setting for First Level Analysis ---%
FirstLevel=1;
if (FirstLevel==1)
    for iSubjID = 1:length(subjects)

        % are extracted from SPM file
        % load event onsets (each subject should have a text file with the BIAC paradigm files arranged according to his/her run order)
        % onsets = textread(['C:\Users\user\Documents\Reyhaneh\EEG+fMRI\EO\MRI_data\Events\' subjects{subjID} '\' analysis '\' subjects{subjID} '_onsets.txt']);
        SubjID=subjects{iSubjID};
        SPM=load([FirstLevelFolder,SubjID,'\3cond_distr2\SPM.mat']);

        % vols = {};
        % are extracted from SPM file
        % run_ind = [1 2 3 4 5 6 7 8]; % we use this to model only the non-SST runs

        nRuns=length(SPM.SPM.Sess);
        disp([SubjID, ': ', num2str(nRuns)])
        for iSess =1:nRuns
            disp(['Session: ',num2str(iSess)])

            bnum=1;
            clear matlabbatch
            matlabbatch{bnum}.spm.stats.fmri_spec.dir = {[rootdir '\' SubjID '\' analysis '\run0' num2str(iSess)]};
            matlabbatch{bnum}.spm.stats.fmri_spec.timing.units = 'scans';
            matlabbatch{bnum}.spm.stats.fmri_spec.timing.RT = 2;
            matlabbatch{bnum}.spm.stats.fmri_spec.timing.fmri_t = 28; % SPM default = 16
            matlabbatch{bnum}.spm.stats.fmri_spec.timing.fmri_t0 = 28; % SPM default = 1

            % are extracted from SPM file
            %         folder = [FuncFolder '\' subjects{subjID} '\run0' num2str(run_ind(r))];
            %         f = dir([folder '\aV*.img']);
            %         vols = {f.name};
            %         nVols = size(vols,2);
            %
            %         for v = 1:nVols
            %             vols{v} = [folder '\' sprintf('aV%04d',v) '.img'];
            %         end
            %
            %         vols = vols';
            % Functional Image:
            FirstSessImgRowID=SPM.SPM.Sess(iSess).row(1);
            LastSessImgRowID=SPM.SPM.Sess(iSess).row(end);

            FirstImg=SPM.SPM.xY.P(FirstSessImgRowID,:);
            RunID=['run',num2str(iSess,'%02.f')];
            BaseEnd=strfind(FirstImg,'Func')+4;
            ImageStart=strfind(FirstImg,'swaV')-1;
            SubjRunFolder=FirstImg(BaseEnd+1:ImageStart);
            FuncName={};
            cntfName=1;
            for iImage=FirstSessImgRowID:LastSessImgRowID

                ImgID=SPM.SPM.xY.P(iImage,ImageStart+3:end);
                MyImage=[FuncFolder,SubjRunFolder,'r',ImgID];%use realigned images (r)
                FuncName{cntfName}=MyImage;
                if ((iImage==FirstSessImgRowID) || (iImage==LastSessImgRowID))
                    if ~exist(MyImage)
                        disp('********************************************')
                        disp('  Error! Func file does not exist!!')
                    end
                end

                cntfName=cntfName+1;
            end
            FuncName=FuncName';
            disp(['   Total functional file: ',num2str(cntfName-1)])
            matlabbatch{bnum}.spm.stats.fmri_spec.sess(1).scans = FuncName;

            R = textread([FuncFolder,SubjRunFolder,'rp_aV0001.txt']); % motion parameters for each run
            %R = load([FuncFolder, SubjRunFolder 'art_regression_outliers_and_movement_aV0001.mat']); % outlier volumes and motion parameters for each run
            ART_r = ([FuncFolder,SubjRunFolder 'art_regression_outliers_raV0001.mat']); % use art output applied on r files

            matlabbatch{bnum}.spm.stats.fmri_spec.sess(1).multi = {''};
            matlabbatch{bnum}.spm.stats.fmri_spec.sess(1).regress = struct('name', {'x_trans' 'y_trans' 'z_trans' 'pitch' 'roll' 'yaw'}, 'val', {R(:,1) R(:,2) R(:,3) R(:,4) R(:,5) R(:,6)});
            %matlabbatch{bnum}.spm.stats.fmri_spec.sess(1).multi_reg = {''};
            %matlabbatch{bnum}.spm.stats.fmri_spec.sess(1).regress = {''};
            matlabbatch{bnum}.spm.stats.fmri_spec.sess(1).multi_reg = {ART_r};
            matlabbatch{bnum}.spm.stats.fmri_spec.sess(1).hpf = 128;

            % Condition Extraction:
            NewCond=struct( 'name',{},'ons',[],'dur',[]);
            cnt=1;
            disp(['  Session_',num2str(iSess)])
            for iCond=1:length(CondTypeName)
                if strcmpi(SPM.SPM.Sess(iSess).U(iCond).name{1},CondTypeName{iCond})
                    for iConNum=1:length(SPM.SPM.Sess(iSess).U(iCond).ons)
                        matlabbatch{bnum}.spm.stats.fmri_spec.sess(1).cond(cnt).name=[SPM.SPM.Sess(iSess).U(iCond).name{1},'_',num2str(iConNum)];
                        matlabbatch{bnum}.spm.stats.fmri_spec.sess(1).cond(cnt).onset=SPM.SPM.Sess(iSess).U(iCond).ons(iConNum);
                        matlabbatch{bnum}.spm.stats.fmri_spec.sess(1).cond(cnt).duration=SPM.SPM.Sess(iSess).U(iCond).dur(iConNum);
                        matlabbatch{bnum}.spm.stats.fmri_spec.sess(1).cond(cnt).tmod = 0;
                        matlabbatch{bnum}.spm.stats.fmri_spec.sess(1).cond(cnt).pmod = struct('name', {}, 'param', {}, 'poly', {});
                        matlabbatch{bnum}.spm.stats.fmri_spec.sess(1).cond(cnt).orth = 1;
                        cnt=cnt+1;
                    end
                    disp(['     ',SPM.SPM.Sess(iSess).U(iCond).name{1},': ',num2str(length(SPM.SPM.Sess(iSess).U(iCond).ons))])
                else
                    disp('********************************************')
                    disp(['      Warning!! Mismatch condition name, ',SPM.SPM.Sess(iSess).U(iCond).name{1},' vs. ',CondTypeName{iCond}])
                end
            end
            TotalCond=cnt-1;
            disp(['   Total: ',num2str(TotalCond)])

            matlabbatch{bnum}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{bnum}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0]; % temp derivative is ON
            matlabbatch{bnum}.spm.stats.fmri_spec.volt = 1;
            matlabbatch{bnum}.spm.stats.fmri_spec.global = 'None';
            matlabbatch{bnum}.spm.stats.fmri_spec.mthresh = 0.8; % SPM default = 0.8
            %matlabbatch{bnum}.spm.stats.fmri_spec.mask = {'Z:\Common\fmri\spm12\tpm\rmask_ICV_4x4x4.nii'};
            matlabbatch{bnum}.spm.stats.fmri_spec.mask = {''};
            matlabbatch{bnum}.spm.stats.fmri_spec.cvi = 'AR(1)';

            bnum = bnum + 1;

            %--- MODEL ESTIMATION ---%
            matlabbatch{bnum}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{bnum-1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat')); % need to specify which bnum this is dependent on
            matlabbatch{bnum}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{bnum}.spm.stats.fmri_est.method.Classical = 1;

            bnum = bnum + 1;

            SPMFile= [matlabbatch{1}.spm.stats.fmri_spec.dir{1},'\SPM.mat'];
            if ~exist(SPMFile)
                spm_jobman('run',matlabbatch);
            end

            % contrast 
            clear matlabbatch
            matlabbatch{1}.spm.stats.con.spmmat = {SPMFile};
            TotalBeta=TotalCond; %3:xyz, 3: raw,pithc,yaw+ 1 constant
            ContrastBase=zeros(1,TotalBeta);
            cnt=1;
            %disp(['  Session_',num2str(iSess)])
            for iCond=1:length(CondTypeName)
                if strcmpi(SPM.SPM.Sess(iSess).U(iCond).name{1},CondTypeName{iCond})
                    for iConNum=1:length(SPM.SPM.Sess(iSess).U(iCond).ons)
                        matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name=[SPM.SPM.Sess(iSess).U(iCond).name{1},'_',num2str(iConNum)];
                        MyWeights=ContrastBase;
                        MyWeights(cnt)=1;
                        matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.weights = MyWeights;
                        matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.sessrep = 'none';
                        cnt=cnt+1;
                    end
%                     disp(['     ',SPM.SPM.Sess(iSess).U(iCond).name{1},': ',num2str(length(SPM.SPM.Sess(iSess).U(iCond).ons))])
                else
                    disp('********************************************')
                    disp(['      Warning!! Mismatch condition name, ',SPM.SPM.Sess(iSess).U(iCond).name{1},' vs. ',CondTypeName{iCond}])
                end
            end   
          
     
            matlabbatch{1}.spm.stats.con.delete = 0;
            spm_jobman('run',matlabbatch)
        end
    end
end

