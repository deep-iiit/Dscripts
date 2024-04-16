% Script README
% SPM based batch processing of fMRI datasets using GLM model
% Author : Deepak Sharma
% Date : October, 2019
% Rename folders names to subject ids 
close all; clear ; clc; % Close windows, Delete all variables and clear screen
%spm fmri % Open SPM FMRI gui for visualization
spm('defaults', 'FMRI'); % For fMRI analysis default variable setting 
delf = [1:9]; % Number of files to delete
max_limit = 1.5;  % Set maximum allowed motion value for realignment
fold_pref = 'sub'; % Folder prefix example: sub001, sub002 ...
cur_d = dir('.'); % PWD
isdir_vec = [cur_d.isdir];
fol_names = {cur_d(isdir_vec).name}';
fol_names(1:2) = [];
total_folders = numel(fol_names);
sid_to_name_mat = cell(total_folders,2);
%%
for i=1:total_folders % Rename folders (Makes them anonymous)
    sid_to_name_mat(i,2) = fol_names(i,1);
    f = pref_set(i,fold_pref,'');
    movefile(fol_names{i}, f)
    sid_to_name_mat(i, 1) = cellstr(f);
end
% Write information to subject_id_to_name_mapping.csv file
map_table = cell2table(sid_to_name_mat);
writetable(map_table, 'subject_id_to_name_mapping.csv')
%% Convert DICOMs to NIFTI files
for i = 1:total_folders
    f = pref_set(i, fold_pref,'');
    cd(f);
    d = dir;
    d = d(~ismember({d.name},{'.','..'}));
    % Select only extensionless files
    matlabbatch{i}.spm.util.import.dicom.data = strcat(f,'/',{d.name}');
    matlabbatch{i}.spm.util.import.dicom.root = 'flat';
    matlabbatch{i}.spm.util.import.dicom.outdir = {pwd};
    matlabbatch{i}.spm.util.import.dicom.protfilter = '.*';
    matlabbatch{i}.spm.util.import.dicom.convopts.format = 'nii';
    matlabbatch{i}.spm.util.import.dicom.convopts.meta = 0;
    matlabbatch{i}.spm.util.import.dicom.convopts.icedims = 0;
    cd ..
end
parfor i=1:total_folders
    spm_jobman('run',matlabbatch(i))
end
save('DICOM_import_batch.mat','matlabbatch')
clear matlabbatch
%%
% Rename structural and functional files
for i=1:total_folders
    f = pref_set(i, fold_pref,'');
    cd(f)
    dirData = dir('f*.nii');
    fileNames = {dirData.name};
    for iFile = 1:numel(fileNames)
      newName = sprintf('fun%03d.nii',iFile);
      movefile(fileNames{iFile},newName);
    end
    dirData = dir('s*.nii');
    fileNames = {dirData.name};
    for iFile = 1:numel(fileNames)
      newName = sprintf('str.nii');
      movefile(fileNames{iFile},newName);
    end
    cd ..
end
%% Remove first few files
for i = 1:total_folders
    f = pref_set(i,fold_pref,'');
    cd(f)
    for j = delf
        f = pref_set(j, 'fun', '.nii');
        delete(f)
    end
    cd ..
end
%% Slice timing correction
for i=1:total_folders
    f = pref_set(i, fold_pref,'');
    cd(f)
    dirData = dir('f*.nii');
    fileNames = {dirData.name}';
    matlabbatch{i}.spm.temporal.st.scans = {strcat(f,'/',fileNames)}';
    matlabbatch{i}.spm.temporal.st.nslices = 36;
    matlabbatch{i}.spm.temporal.st.tr = 3;
    matlabbatch{i}.spm.temporal.st.ta = 2.9166;
    matlabbatch{i}.spm.temporal.st.so = [2:2:36 1:2:36];
    matlabbatch{i}.spm.temporal.st.refslice = 18;
    matlabbatch{i}.spm.temporal.st.prefix = 'a';
    cd ..
end
for i=1:total_folders
    spm_jobman('run',matlabbatch(i))
end
save('Slice_time_correction_batch.mat','matlabbatch')
clear matlabbatch
%% Realignment of functional files
for i = 1:total_folders
    f = pref_set(i, fold_pref,'');
    cd(f)
    dirData = dir('af*.nii');
    fileNames = {dirData.name}';
    matlabbatch{i}.spm.spatial.realign.estwrite.data = {strcat(f,'\',fileNames)}';
    matlabbatch{i}.spm.spatial.realign.estwrite.eoptions.quality = 1;
    matlabbatch{i}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{i}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{i}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    matlabbatch{i}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{i}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{i}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{i}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    matlabbatch{i}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{i}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{i}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{i}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    cd ..
end
pwd
for i=1:total_folders
    f = pref_set(i, fold_pref,'');
    spm_jobman('run',matlabbatch(i))
end
save('Realignment_batch.mat','matlabbatch')
clear matlabbatch
%% Convert 6 to 24 Friston regressors
for i=1:total_folders
    fprintf('\tConverting 6 to 24 Friston regressors for subject no.%d\n',i);
    f = pref_set(i,fold_pref,'');
    cd(f);
    filename = dir('rp_afun*.txt');
    in = readmatrix(filename.name);
    format('long','e')
    rows = size(in, 1);
    out = zeros(rows, 24);
    for r = 1:rows
        if r ~= 1
            out(r, 1)  = in(r, 1) ;
            out(r, 2)  = in(r, 2) ;
            out(r, 3)  = in(r, 3) ;
            out(r, 4)  = in(r, 4) ;
            out(r, 5)  = in(r, 5) ;
            out(r, 6)  = in(r, 6) ;
            out(r, 7)  = in(r, 1) ^ 2 ;
            out(r, 8)  = in(r, 2) ^ 2 ;
            out(r, 9)  = in(r, 3) ^ 2 ;
            out(r, 10) = in(r, 4) ^ 2 ;
            out(r, 11) = in(r, 5) ^ 2 ;
            out(r, 12) = in(r, 6) ^ 2 ;
            out(r, 13) = in(r, 1) - in(r - 1, 1) ;
            out(r, 14) = in(r, 2) - in(r - 1, 2) ;
            out(r, 15) = in(r, 3) - in(r - 1, 3) ;
            out(r, 16) = in(r, 4) - in(r - 1, 4) ;
            out(r, 17) = in(r, 5) - in(r - 1, 5) ;
            out(r, 18) = in(r, 6) - in(r - 1, 6) ;
            out(r, 19) = out(r, 13) ^ 2 ;
            out(r, 20) = out(r, 14) ^ 2 ;
            out(r, 21) = out(r, 15) ^ 2 ;
            out(r, 22) = out(r, 16) ^ 2 ;
            out(r, 23) = out(r, 17) ^ 2 ;
            out(r, 24) = out(r, 18) ^ 2 ;
        end   
    end
    out(1,:) = 0;
    save regressors.txt out -ascii -tabs
    present_path = pwd;
    splits = strsplit(present_path, '/');
    subj_folder = splits{end};
    reg_filename = strcat(subj_folder, '_regressors.txt');
    movefile('regressors.txt', reg_filename);
    cd ..
end
%% Summarize realignment results
disp('Testing quality of data ...');
format long
quality = zeros(total_folders,14);
for i=1:total_folders
    fprintf('\tTesting for subject no.%d\n',i);
    f = pref_set(i,fold_pref,'');
    cd(f);
    filename = dir('rp_afun*.txt');
    fileID = fopen(filename.name,'r');
    C = textscan(fileID,'%f %f %f %f %f %f','Delimiter',' ','MultipleDelimsAsOne',1);
    fclose(fileID);
    a = cell2mat(C);
    quality(i,1) = i;
    quality(i,2) = max(abs(a(:,1)));
    quality(i,3) = max(abs(a(:,2)));
    quality(i,4) = max(abs(a(:,3)));
    quality(i,5) = max(abs(a(:,4)));
    quality(i,6) = max(abs(a(:,5)));
    quality(i,7) = max(abs(a(:,6)));
    if(quality(i,2) > max_limit || quality(i,3) > max_limit || quality(i,4) > max_limit || quality(i,5) > max_limit || quality(i,6) > max_limit || quality(i,7) > max_limit)
        quality(i,8) = 1;
    end
    quality(i,9)  = std(a(:,1));
    quality(i,10) = std(a(:,2));
    quality(i,11) = std(a(:,3));
    quality(i,12) = std(a(:,4));
    quality(i,13) = std(a(:,5));
    quality(i,14) = std(a(:,6));
    cd ..
end
fid = fopen('quality.csv', 'w+');   % Write information to quality.csv file
for i=1:size(quality, 1)
    fprintf(fid, '%f ,', quality(i,:));
    fprintf(fid, '\n');
end
fclose(fid);
fprintf('Result file path: %s',strcat(pwd,'/quality.csv'));
fprintf('\nTesting done ...\n');
%% Coregistration
clear matlabbatch
for i = 1:total_folders
    f = pref_set(i, fold_pref,'');
    matlabbatch{i}.spm.spatial.coreg.estimate.ref = {['./' f '/meanafun010.nii,1']};
    matlabbatch{i}.spm.spatial.coreg.estimate.source = {['./' f '/str.nii,1']};
    matlabbatch{i}.spm.spatial.coreg.estimate.other = {{}}';
    matlabbatch{i}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{i}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{i}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{i}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
end
for i=1:total_folders
    spm_jobman('run',matlabbatch(i))
end
save('Coregistration_batch.mat','matlabbatch')
clear matlabbatch
%% Segmentation
clear matlabbatch
for i=1:total_folders
    f = pref_set(i, fold_pref,'');
    matlabbatch{i}.spm.spatial.preproc.channel.vols = {['./' f '/str.nii,1']};
    matlabbatch{i}.spm.spatial.preproc.channel.biasreg = 1.000000000000000e-03;
    matlabbatch{i}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{i}.spm.spatial.preproc.channel.write = [0,1];
    matlabbatch{i}.spm.spatial.preproc.tissue(1).tpm = {'/home/deep/matlab/spm12/tpm/TPM.nii,1'};
    matlabbatch{i}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{i}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{i}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{i}.spm.spatial.preproc.tissue(2).tpm = {'/home/deep/matlab/spm12/tpm/TPM.nii,2'};
    matlabbatch{i}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{i}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{i}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{i}.spm.spatial.preproc.tissue(3).tpm = {'/home/deep/matlab/spm12/tpm/TPM.nii,3'};
    matlabbatch{i}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{i}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{i}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{i}.spm.spatial.preproc.tissue(4).tpm = {'/home/deep/matlab/spm12/tpm/TPM.nii,4'};
    matlabbatch{i}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{i}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{i}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{i}.spm.spatial.preproc.tissue(5).tpm = {'/home/deep/matlab/spm12/tpm/TPM.nii,4'};
    matlabbatch{i}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{i}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{i}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{i}.spm.spatial.preproc.tissue(6).tpm = {'/home/deep/matlab/spm12/tpm/TPM.nii,6'};
    matlabbatch{i}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{i}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{i}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{i}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{i}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{i}.spm.spatial.preproc.warp.reg = [0,1.000000000000000e-03,0.500000000000000,0.050000000000000,0.200000000000000];
    matlabbatch{i}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{i}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{i}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{i}.spm.spatial.preproc.warp.write = [0,1];
end
for i=1:total_folders
    spm_jobman('run',matlabbatch(i))
end
save('Segmentation_batch.mat','matlabbatch')
clear matlabbatch
%% Normalization
clc; clear matlabbatch
for i=1:total_folders
    f = pref_set(i, fold_pref,'');
    cd(f)
    dirData = dir('raf*.nii');
    fileNames = {dirData.name}';
    matlabbatch{i}.spm.spatial.normalise.write.subj.def = {['./' f '/y_str.nii']};
    matlabbatch{i}.spm.spatial.normalise.write.subj.resample = strcat(f,'/',fileNames);
    matlabbatch{i}.spm.spatial.normalise.write.woptions.bb = [-78,-112,-70;78,76,85];
    matlabbatch{i}.spm.spatial.normalise.write.woptions.vox = [3,3,3];
    matlabbatch{i}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{i}.spm.spatial.normalise.write.woptions.prefix = 'w';
    cd ..
end
for i=1:total_folders
    spm_jobman('run',matlabbatch(i))
end
save('Normalization_batch.mat','matlabbatch')
clear matlabbatch
disp('Normalization complete')
%% Smoothing
clc; clear matlabbatch
for i=1:total_folders
    f = pref_set(i, fold_pref,'');
    cd(f)
    dirData = dir('wraf*.nii');
    fileNames = {dirData.name}';
    matlabbatch{i}.spm.spatial.smooth.data = strcat(f,'/',fileNames);
    matlabbatch{i}.spm.spatial.smooth.fwhm =[6,6,6];
    matlabbatch{i}.spm.spatial.smooth.dtype = 0;
    matlabbatch{i}.spm.spatial.smooth.im = 0;
    matlabbatch{i}.spm.spatial.smooth.prefix = 's';
    cd ..
end
for i=1:total_folders
    spm_jobman('run',matlabbatch(i))
end
save('Smoothing_batch.mat','matlabbatch')
clear matlabbatch
disp('Smoothing complete')
%%
% First level specification
for i=1:total_folders
    f = pref_set(i, fold_pref,'');
    cd(f)
    dirData = dir('swraf*.nii');
    fileNames = {dirData.name}';
    regData = dir('rp_afun*.txt');
    reg_files = {regData.name}';
    matlabbatch{i}.spm.stats.fmri_spec.dir = {f};
    matlabbatch{i}.spm.stats.fmri_spec.timing.units = 'scans';
    matlabbatch{i}.spm.stats.fmri_spec.timing.RT = 3;
    matlabbatch{i}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{i}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    matlabbatch{i}.spm.stats.fmri_spec.sess.scans = strcat(f,'\',fileNames);
    matlabbatch{i}.spm.stats.fmri_spec.sess.cond(1).name = 'DL_baseline';
    matlabbatch{i}.spm.stats.fmri_spec.sess.cond(1).onset = [19 47 75];
    matlabbatch{i}.spm.stats.fmri_spec.sess.cond(1).duration = [9 9 9];
    matlabbatch{i}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
    matlabbatch{i}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{i}.spm.stats.fmri_spec.sess.cond(1).orth = 1;
    matlabbatch{i}.spm.stats.fmri_spec.sess.cond(2).name = 'DL_active';
    matlabbatch{i}.spm.stats.fmri_spec.sess.cond(2).onset = [0 28 56];
    matlabbatch{i}.spm.stats.fmri_spec.sess.cond(2).duration = [19 19 19]';
    matlabbatch{i}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
    matlabbatch{i}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{i}.spm.stats.fmri_spec.sess.cond(2).orth = 1;
    matlabbatch{i}.spm.stats.fmri_spec.sess.cond(3).name = 'NV_baseline';
    matlabbatch{i}.spm.stats.fmri_spec.sess.cond(3).onset = [103 131 159];
    matlabbatch{i}.spm.stats.fmri_spec.sess.cond(3).duration = [9 9 9];
    matlabbatch{i}.spm.stats.fmri_spec.sess.cond(3).tmod = 0;
    matlabbatch{i}.spm.stats.fmri_spec.sess.cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{i}.spm.stats.fmri_spec.sess.cond(3).orth = 1;
    matlabbatch{i}.spm.stats.fmri_spec.sess.cond(4).name = 'NV_active';
    matlabbatch{i}.spm.stats.fmri_spec.sess.cond(4).onset = [84 112 140];
    matlabbatch{i}.spm.stats.fmri_spec.sess.cond(4).duration = [19 19 19]';
    matlabbatch{i}.spm.stats.fmri_spec.sess.cond(4).tmod = 0;
    matlabbatch{i}.spm.stats.fmri_spec.sess.cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{i}.spm.stats.fmri_spec.sess.cond(4).orth = 1;
    matlabbatch{i}.spm.stats.fmri_spec.sess.multi = {''};
    matlabbatch{i}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
    matlabbatch{i}.spm.stats.fmri_spec.sess.multi_reg = strcat(f,'\',reg_files);
    matlabbatch{i}.spm.stats.fmri_spec.sess.hpf = 128;
    matlabbatch{i}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{i}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{i}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{i}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{i}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{i}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{i}.spm.stats.fmri_spec.cvi = 'AR(1)';
    cd ..
end
for i=1:total_folders
    spm_jobman('run',matlabbatch(i))
end
save('First_level_specification_batch.mat','matlabbatch')
clear matlabbatch
% First level estimation
for i=1:total_folders
    f = pref_set(i, fold_pref,'');
    cd(f)
    dirData = dir('SPM.mat');
    fileNames = {dirData.name}';
    matlabbatch{i}.spm.stats.fmri_est.spmmat = strcat(f,'\',fileNames);
    matlabbatch{i}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{i}.spm.stats.fmri_est.method.Classical = 1;
    cd ..
end
for i=1:total_folders
    spm_jobman('run',matlabbatch(i))
end
save('First_level_Estimation_batch.mat','matlabbatch')
clear matlabbatch
% First level contrasts
for i=1:total_folders
    f = pref_set(i, fold_pref,'');
    cd(f)
    dirData = dir('SPM.mat');
    fileNames = {dirData.name}';
    matlabbatch{i}.spm.stats.con.spmmat = strcat(f,'\',fileNames);
    matlabbatch{i}.spm.stats.con.consess{1}.tcon.name = 'DLbaseline-DLactive';
    matlabbatch{i}.spm.stats.con.consess{1}.tcon.weights = [1 -1 0 0];
    matlabbatch{i}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.consess{2}.tcon.name = 'DLactive-DLbaseline';
    matlabbatch{i}.spm.stats.con.consess{2}.tcon.weights = [-1 1 0 0];
    matlabbatch{i}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.consess{3}.tcon.name = 'NVbaseline-NVactive';
    matlabbatch{i}.spm.stats.con.consess{3}.tcon.weights = [0 0 1 -1];
    matlabbatch{i}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.consess{4}.tcon.name = 'NVactive-NVbaseline';
    matlabbatch{i}.spm.stats.con.consess{4}.tcon.weights = [0 0 -1 1];
    matlabbatch{i}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.consess{5}.tcon.name = 'NVactive-DLactive';
    matlabbatch{i}.spm.stats.con.consess{5}.tcon.weights = [0 -1 0 1];
    matlabbatch{i}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.consess{6}.tcon.name = 'DLactive-NVactive';
    matlabbatch{i}.spm.stats.con.consess{6}.tcon.weights = [0 1 0 -1];
    matlabbatch{i}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.consess{7}.tcon.name = 'NVbaseline-DLbaseline';
    matlabbatch{i}.spm.stats.con.consess{7}.tcon.weights = [-1 0 1 0];
    matlabbatch{i}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.consess{8}.tcon.name = 'DLbaseline-NVbaseline';
    matlabbatch{i}.spm.stats.con.consess{8}.tcon.weights = [1 0 -1 0];
    matlabbatch{i}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
    matlabbatch{i}.spm.stats.con.delete = 1;
    cd ..
end
for i=1:total_folders
    spm_jobman('run',matlabbatch(i))
end
save('Contrast_batch.mat','matlabbatch')
ncons = length(matlabbatch{1,1}.spm.stats.con.consess);
clear matlabbatch
% First level contrast results
for j = 1:ncons
    for i=1:total_folders
        f = pref_set(i, fold_pref,'');
        cd(f)
        dirData = dir('SPM.mat');
        fileNames = {dirData.name}';
        matlabbatch{i}.spm.stats.results.spmmat = strcat(pwd,'\',fileNames);
        matlabbatch{i}.spm.stats.results.conspec.titlestr = '';
        matlabbatch{i}.spm.stats.results.conspec.contrasts = j;
        matlabbatch{i}.spm.stats.results.conspec.threshdesc = 'FWE';
        matlabbatch{i}.spm.stats.results.conspec.thresh = 0.05;
        matlabbatch{i}.spm.stats.results.conspec.extent = 0;
        matlabbatch{i}.spm.stats.results.conspec.conjunction = 1;
        matlabbatch{i}.spm.stats.results.conspec.mask.none = 1;
        matlabbatch{i}.spm.stats.results.units = 1;
        matlabbatch{i}.spm.stats.results.export{1}.png = true;
        cd ..
    end
    for i=1:total_folders
        spm_jobman('run',matlabbatch(i))
    end
    save('Results_batch.mat','matlabbatch')
    clear matlabbatch
    cd ..
end
%% Second level analysis - One sample t-test
mkdir group_analysis_one_sample_t_test
clc
wd = pwd;
cons = {};
% cfile = 'con_0001.nii'; conname = 'DLbaseline-DLactive';
% cfile = 'con_0002.nii'; conname = 'DLactive-DLbaseline';
% cfile = 'con_0003.nii'; conname = 'NVbaseline-NVactive';
% cfile = 'con_0004.nii'; conname = 'NVactive-NVbaseline';
% cfile = 'con_0005.nii'; conname = 'NVactive-DLactive';
% cfile = 'con_0006.nii'; conname = 'DLactive-NVactive';
% cfile = 'con_0007.nii'; conname = 'NVbaseline-DLbaseline';
cfile = 'con_0008.nii'; conname = 'DLbaseline-NVbaseline';
for i=1:total_folders
    f = pref_set(i, fold_pref,'');
    cd(f)
    dirData = dir(cfile);
    confile = {dirData.name}';
    cons{i} = (strcat(wd,'\',f,'\',confile,',1'));
    cd ..
end
cons = cons';
for i=1:total_folders
    cons{i} = cons{i, 1}{1, 1};
end
matlabbatch{1}.spm.stats.factorial_design.dir = {strcat(pwd,'\group_analysis_one_sample_t_test')};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = cons;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
spm_jobman('run',matlabbatch)
save('group_analysis_one_sample_t_test_design.mat','matlabbatch')
clear matlabbatch
% Estimate for group_analysis_one_sample_t_test
cd group_analysis_one_sample_t_test
dirData = dir('SPM.mat');
SPMfile = {dirData.name}';
matlabbatch{1}.spm.stats.fmri_est.spmmat = SPMfile;
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run',matlabbatch)
cd ..
save('group_analysis_one_sample_t_test_design_estimation.mat','matlabbatch')
clear matlabbatch
% Group_level_one_sample_contrast
cd group_analysis_one_sample_t_test
dirData = dir('SPM.mat');
SPMfile = {dirData.name}';
matlabbatch{1}.spm.stats.con.spmmat = strcat(pwd,'\',SPMfile);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = conname;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;
spm_jobman('run',matlabbatch)
cd ..
save('Group_level_one_sample_contrast.mat','matlabbatch')
clear matlabbatch
% Group_level_one_sample_result
cd group_analysis_one_sample_t_test
dirData = dir('SPM.mat');
SPMfile = {dirData.name}';
matlabbatch{1}.spm.stats.results.spmmat = strcat(pwd,'\',SPMfile);
matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
matlabbatch{1}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{1}.spm.stats.results.conspec.extent = 0;
matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.export{1}.png = true;
spm_jobman('run',matlabbatch)
cd ..
save('Group_level_one_sample_result.mat','matlabbatch')
clear matlabbatch
%% Group analysis - Two sample t-test
mkdir group_analysis_two_sample_t_test
group1 = [12,20,2,10,18,4,13,6,19,5,15];
group2 = [14,16,1,7,9,11,17,21,3,8];
scans1 = {};
scans2 = {};
j = 1; k = 1;
wd = pwd;
% cfile = 'con_0001.nii'; conname = 'DLbaseline-DLactive';
% cfile = 'con_0002.nii'; conname = 'DLactive-DLbaseline';
% cfile = 'con_0003.nii'; conname = 'NVbaseline-NVactive';
% cfile = 'con_0004.nii'; conname = 'NVactive-NVbaseline';
cfile = 'con_0005.nii'; conname = 'NVactive-DLactive';
% cfile = 'con_0006.nii'; conname = 'DLactive-NVactive';
% cfile = 'con_0007.nii'; conname = 'NVbaseline-DLbaseline';
% cfile = 'con_0008.nii'; conname = 'DLbaseline-NVbaseline';
for i=1:total_folders
    f = pref_set(i, fold_pref,'');
    cd(f)
    if(ismember(i,group1))
        dirData = dir(cfile);
        scanfile = {dirData.name}';
        scans1{j} = (strcat(wd,'\',f,'\',scanfile,',1'));
        j = j + 1;
    elseif(ismember(i,group2))
        dirData = dir(cfile);
        scanfile = {dirData.name}';
        scans2{k} = (strcat(wd,'\',f,'\',scanfile,',1'));
        k = k + 1;
    end
    cd ..
end
scans1 = scans1';
scans2 = scans2';
for j = 1:length(scans1)
    scans1{j} = scans1{j,1}{1,1};
end
for j = 1:length(scans2)
    scans2{j} = scans2{j,1}{1,1};
end
matlabbatch{1}.spm.stats.factorial_design.dir = {strcat(pwd,'\group_analysis_two_sample_t_test')};
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = scans1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = scans2;
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
cd group_analysis_two_sample_t_test
spm_jobman('run',matlabbatch)
cd ..
save('Group_level_two_sample_design.mat','matlabbatch')
clear matlabbatch
% Two_sample_t_test_estimate
cd group_analysis_two_sample_t_test
dirData = dir('SPM.mat');
SPMfile = {dirData.name}';
matlabbatch{1}.spm.stats.fmri_est.spmmat = strcat(pwd,'\',SPMfile);
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run',matlabbatch)
cd ..
save('Group_level_two_sample_estimate.mat','matlabbatch')
clear matlabbatch
% Two_sample_t_test_contrast
cd group_analysis_two_sample_t_test
dirData = dir('SPM.mat');
SPMfile = {dirData.name}';
matlabbatch{1}.spm.stats.con.spmmat = strcat(pwd,'\',SPMfile);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = conname;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = conname;
matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;
spm_jobman('run',matlabbatch)
cd ..
save('Group_level_two_sample_contrast.mat','matlabbatch')
clear matlabbatch
% Two_sample_t_test_result
cd group_analysis_two_sample_t_test
dirData = dir('SPM.mat');
SPMfile = {dirData.name}';
matlabbatch{1}.spm.stats.results.spmmat = strcat(pwd,'\',SPMfile);
matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
matlabbatch{1}.spm.stats.results.conspec.contrasts = Inf;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{1}.spm.stats.results.conspec.extent = 0;
matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.export{1}.png = true;
spm_jobman('run',matlabbatch)
cd ..
save('Group_level_two_sample_t_test_result.mat','matlabbatch')
clear matlabbatch






