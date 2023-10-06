% dg_fcmri_wl_wholeGroup
clear all; close all;
cd('/home/uni10/nmri/projects/dgarnica/MRI_EEG_PROSPECTIVE');
% run('dg_open_Conn_19b.m'); % LOADS SPM12 AND CONN
addpath(genpath('updateThisDirectory/spm12_7487'));
addpath(genpath('updateThisDirectory/conn_19b'));
conn_folder = '/home/uni10/nmri/tools/conn/conn_19b';

%% FIND AND ADD functional/structural files
% note: this will look for all data in these folders, irrespective of the specific download subsets entered as command-line arguments

data_dir = '/home/uni10/nmri/projects/dgarnica/MRI_EEG_PROSPECTIVE';
cwd = '/home/uni10/nmri/projects/dgarnica/MRI_EEG_PROSPECTIVE/fcFMRI_analysis/wl_regions';
addpath(fullfile(data_dir, 'fcFMRI_analysis'));

patients_list = {'pseudonyms'};
controls_list = {'pseudonyms'};

% patient's functional files (UNPROCESSED)
pf_files = [];
for i = 1:length(patients_list)
    pf_files{i} = fullfile(data_dir, 'ROLANDIC', patients_list{i}, 'func', '--dataset--', [patients_list{i}, '.EPI_fMRI.language.PA.nii']);
end
pf_files = pf_files';

% control's functional files (UNPROCESSED)
cf_files = [];
for i = 1:length(controls_list)
    cf_files{i} = fullfile(data_dir, 'CONTROLS', controls_list{i}, 'func', '--dataset--', [controls_list{i}, '.EPI_fMRI.language.PA.nii']);
end
cf_files = cf_files';

% patient's structural files
ps_files = [];
for i = 1:length(patients_list)
    ps_files{i} = fullfile(data_dir, 'ROLANDIC', patients_list{i}, 'anat', [patients_list{i}, '.hrT1.nii']);
end
ps_files = ps_files';

% control's structural files
cs_files = [];
for i = 1:length(controls_list)
    cs_files{i} = fullfile(data_dir, 'CONTROLS', controls_list{i}, 'anat', [controls_list{i}, '.hrT1.nii']);
end
cs_files = cs_files';

nsubjects = length(patients_list) + length(controls_list);
functional_file = cat(1,pf_files, cf_files);
structural_file = cat(1,ps_files, cs_files);

if rem(length(functional_file), nsubjects)
    error('mismatch number of functional files %n', length(functional_file));
end

if rem(length(structural_file), nsubjects)
    error('mismatch number of anatomical files %n', length(functional_file));
end

nsessions = length(functional_file) / nsubjects;
functional_file = reshape(functional_file, [nsessions, nsubjects]);
structural_file = {structural_file{1:nsubjects}};

disp([num2str(size(functional_file, 1)), ' sessions']);
disp([num2str(size(functional_file, 2)), ' subjects']);
TR = 1.4; % Repetition time

%% CONN-SPECIFIC SECTION: RUNS PREPROCESSING/SETUP/DENOISING/ANALYSIS STEPS
% Prepares batch structure
clear batch;
batch.filename=fullfile(cwd,'wl_regions.mat');                   % New conn_*.mat experiment name

% SETUP & PREPROCESSING step (using default values for most parameters, see help conn_batch to define non-default values)
% CONN Setup                                             % Default options (uses all ROIs in conn/rois/ directory); see conn_batch for additional options
% CONN Setup.preprocessing                               (realignment/coregistration/segmentation/normalization/smoothing)
batch.Setup.isnew=1;
batch.Setup.nsubjects=nsubjects;
batch.Setup.RT=TR;                                       % TR (seconds)
batch.Setup.functionals=repmat({{}},[nsubjects,1]);      % Point to functional volumes for each subject/session

for nsub = 1:nsubjects
    for nses = 1:nsessions
        batch.Setup.functionals{nsub}{nses}{1} = functional_file{nses, nsub};
    end
end

batch.Setup.structurals = structural_file; % Point to anatomical volumes for each subject
batch.Setup.conditions = struct;
batch.Setup.conditions.names = {'Animals', 'Fractals'}; % 2 conditions within the language task: Animals (phonological) and Fractals (visual)
nconditions = length(batch.Setup.conditions.names); 

subjects_list = cat(2,patients_list,controls_list);
for nses = 1:nsessions
    for ncond = 1:nconditions
        for nsub = 1:nsubjects
            clear subject;
            % call proper conf_file, to get ONSETS and DURATIONS (Durations are the same Reaction Times, in WL task)
            if nsub <= size(patients_list,2)
                conf_file = fullfile(data_dir, 'ROLANDIC', subjects_list{nsub}, 'conf', 'EPI_fMRI_language_PA_mcf_topped_nii', 'subjectinfo_EPI_fMRI_language_PA_mcf_topped_nii.mat');
                load(conf_file);
            elseif nsub > size(patients_list,2)
                conf_file = fullfile(data_dir, 'CONTROLS', subjects_list{nsub}, 'conf', 'EPI_fMRI_language_PA_mcf_topped_nii', 'subjectinfo_EPI_fMRI_language_PA_mcf_topped_nii.mat');
                load(conf_file);                
            end
            
            % check if onsets and durations must be trimmed (if some are greater than dataset length)
            [subject] = dg_child_taskdata_trim_mod(subject);            
            % assign proper onsets and durations (Animals or Fractals)
            if ncond == 1
                batch.Setup.conditions.onsets{ncond}{nsub}{nses} = [subject.WLtask.Onsets_animals];
                batch.Setup.conditions.durations{ncond}{nsub}{nses} = [subject.WLtask.RT_animals];
            elseif ncond == 2
                batch.Setup.conditions.onsets{ncond}{nsub}{nses} = [subject.WLtask.Onsets_fractals];
                batch.Setup.conditions.durations{ncond}{nsub}{nses} = [subject.WLtask.RT_fractals];
            end
        end
    end
end % session-specific conditions

% SETUP ROIs: STOP HERE AND SELECT REGIONS/NETWORKS TO ANALYZE
% batch.Setup.rois.files = {fullfile(conn_folder, 'rois', 'atlas.nii')};
batch.Setup.rois.files{1} = {fullfile(data_dir, 'fcFMRI_analysis', 'rois_images', 'phonological_uniformity-test_z_FDR_0.01.nii')}; % from https://neurosynth.org/analyses/terms/phonological/
batch.Setup.rois.files{2} = {fullfile(data_dir, 'fcFMRI_analysis', 'rois_images', 'reading_uniformity-test_z_FDR_0.01.nii')}; % from https://neurosynth.org/analyses/terms/reading/
load(fullfile(data_dir, 'fcFMRI_analysis', 'ROIs')); % load the ROIs I saved from conn_folder = '/home/uni10/nmri/tools/conn/conn_19b';
% Language and reading related regions
batch.Setup.rois.names = {'FP r','FP l','IC r','IC l','SFG r','SFG l','MidFG r','MidFG l','IFG tri r','IFG tri l','IFG oper r','IFG oper l','PreCG r','PreCG l','TP r','TP l','aSTG r','aSTG l','pSTG r','pSTG l','aMTG r','aMTG l','pMTG r','pMTG l','toMTG r','toMTG l','aITG r','aITG l','pITG r','pITG l','toITG r','toITG l','PostCG r','PostCG l','aSMG r','aSMG l','pSMG r','pSMG l','AG r','AG l','LG r','LG l','aTFusC r','aTFusC l','pTFusC r','pTFusC l','TOFusC r','TOFusC l','OFusG r','OFusG l','HG r','HG l','PT r','PT l'};
batch.Setup.rois.networks = 0; % Ensure to turn off network-level analysis

batch.Setup.preprocessing.steps='default_mni'; % RUNS DEFAULT MNI PROCESSING, AVAILABLE IN CONN
% batch.Setup.preprocessing.steps = setdiff(batch.Setup.preprocessing.steps, {'functional_realign&unwarp', 'functional_slice-timing'}); % Remove Realignment & Slice-timing correction steps (SINCE DATA ALREADY PREPROCESSED WITH TOPUP+MCFLIRT WAS SELECTED)
batch.Setup.preprocessing.sliceorder='interleaved (Siemens)';
batch.Setup.done=1;
batch.Setup.overwrite='Yes';

% Run one step at a time:
% conn_batch(batch); % runs Preprocessing and Setup steps only
% clear batch;
% batch.filename=fullfile(cwd,'wl_regions.mat');            % Existing conn_*.mat experiment name

%% DENOISING step
% CONN Denoising                                    % Default options (uses White Matter+CSF+realignment+scrubbing+conditions as confound regressors); see conn_batch for additional options
batch.Denoising.filter=[0.01, 0.1];                 % frequency filter (band-pass values, in Hz)
batch.Denoising.done=1;
batch.Denoising.overwrite='Yes';

% Run one step at a time:
% conn_batch(batch); % runs Denoising step only 
% clear batch;
% batch.filename=fullfile(cwd,'wl_regions.mat');            % Existing conn_*.mat experiment name

%% FIRST-LEVEL ANALYSIS (SCRIPTED) AND SECOND-LEVEL ANALYSIS (TO DO IN GUI)
% **First-level Analysis - Default options for (except that uses previously defined ROIs as connectivity sources)
batch.Analysis.done = 1;
batch.Analysis.overwrite = 'Yes';
batch.Analysis.name = 'S2V_R2R'; % Seed-to-Voxel & Weighted ROI-to-ROI
batch.Analysis.measure = 1;     % connectivity measure used, 1='correlation (bivariate)', 2='correlation (semipartial)', 
                                % 3='regression (bivariate)', 4='regression (multivariate)'; [1] 
batch.Analysis.weight = 2;      % within-condition weight, 1 = 'none', 2 = 'hrf', 3 = 'hanning'; [2]
batch.Analysis.modulation = 0;  % temporal modulation, 0 = standard weighted GLM analyses; 1 = gPPI analyses of condition-specific 
                                % temporal modulation factor, or a string for PPI analyses of other temporal modulation factor 
                                % (same for all conditions; valid strings are ROI names and 1st-level covariate names)'; [0]
batch.Analysis.conditions = []; % (for modulation==1 only) list of task condition names to be simultaneously entered in gPPI 
                                % model (leave empty for default 'all existing conditions') []
batch.Analysis.type = 3;        % analysis type, 1 = 'ROI-to-ROI', 2 = 'Seed-to-Voxel', 3 = 'all'; [3]

% Run ALL analyses
conn_batch(batch);

% CONN Display
% launches conn gui to explore results AND RUN 2nd LEVEL ANALYSIS WITH GUI
conn
conn('load', fullfile(cwd, 'wl_regions.mat'));
conn gui_results
