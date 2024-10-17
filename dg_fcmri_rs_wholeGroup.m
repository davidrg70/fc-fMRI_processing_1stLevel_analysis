% dg_fcmri_rs_wholeGroup
clear all; close all; clc;
cd('PROJECT/DIRECTORY');
% run('dg_open_Conn_19b.m'); % LOADS SPM12 AND CONN
addpath(genpath('updateThisDirectory/spm12_7487'));
addpath(genpath('updateThisDirectory/conn_19b'));

%% FIND AND ADD functional/structural files
% note: this will look for all data in these folders, irrespective of the specific download subsets entered as command-line arguments
data_dir = 'PROJECT/DIRECTORY';
cwd = 'RESULTS/DIRECTORY/fcFMRI_analysis/rs_regions';
addpath(fullfile(data_dir, 'fcFMRI_analysis'));

patients_list = {'pseudonyms'};
controls_list = {'pseudonyms'};

% patient's functional files (UNPROCESSED)
pf_files = [];
for i = 1:length(patients_list)
    pf_files{i} = fullfile(data_dir, 'ROLANDIC', patients_list{i}, 'func', '--dataset--', [patients_list{i}, '.EPI_fMRI.rs.PA.nii']);
end
pf_files = pf_files';

% control's functional files (UNPROCESSED)
cf_files = [];
for i = 1:length(controls_list)
    cf_files{i} = fullfile(data_dir, 'CONTROLS', controls_list{i}, 'func', '--dataset--', [controls_list{i}, '.EPI_fMRI.rs.PA.nii']);
end
cf_files = cf_files';

% patient's structural files (UNPROCESSED)
ps_files = [];
for i = 1:length(patients_list)
    ps_files{i} = fullfile(data_dir, 'ROLANDIC', patients_list{i}, 'anat', [patients_list{i}, '.hrT1.nii']);
end
ps_files = ps_files';

% control's structural files (UNPROCESSED)
cs_files = [];
for i = 1:length(controls_list)
    cs_files{i} = fullfile(data_dir, 'CONTROLS', controls_list{i}, 'anat', [controls_list{i}, '.hrT1.nii']);
end
cs_files = cs_files';

NSUBJECTS = length(patients_list) + length(controls_list);
FUNCTIONAL_FILE = cat(1,pf_files, cf_files);
STRUCTURAL_FILE = cat(1,ps_files, cs_files);

if rem(length(FUNCTIONAL_FILE), NSUBJECTS)
    error('mismatch number of functional files %n', length(FUNCTIONAL_FILE));
end

if rem(length(STRUCTURAL_FILE), NSUBJECTS)
    error('mismatch number of anatomical files %n', length(FUNCTIONAL_FILE));
end

nsessions = length(FUNCTIONAL_FILE) / NSUBJECTS;
FUNCTIONAL_FILE = reshape(FUNCTIONAL_FILE, [nsessions, NSUBJECTS]);
STRUCTURAL_FILE = {STRUCTURAL_FILE{1:NSUBJECTS}};

disp([num2str(size(FUNCTIONAL_FILE, 1)), ' sessions']);
disp([num2str(size(FUNCTIONAL_FILE, 2)), ' subjects']);
TR = 1.4; % Repetition time

%% CONN-SPECIFIC SECTION: RUNS PREPROCESSING/SETUP/DENOISING/ANALYSIS STEPS
% Prepares batch structure
clear batch;
batch.filename=fullfile(cwd,'rs_regions.mat');            % New conn_*.mat experiment name

% SETUP & PREPROCESSING step (using default values for most parameters, see help conn_batch to define non-default values)
% CONN Setup                                            % Default options (uses all ROIs in conn/rois/ directory); see conn_batch for additional options
% CONN Setup.preprocessing                               (realignment/coregistration/segmentation/normalization/smoothing)
batch.Setup.isnew=1;
batch.Setup.nsubjects=NSUBJECTS;
batch.Setup.RT=TR;                                        % TR (seconds)
batch.Setup.functionals=repmat({{}},[NSUBJECTS,1]);       % Point to functional volumes for each subject/session

for nsub = 1:NSUBJECTS
    for nses = 1:nsessions
        batch.Setup.functionals{nsub}{nses}{1} = FUNCTIONAL_FILE{nses, nsub};
    end
end

batch.Setup.structurals = STRUCTURAL_FILE; % Point to anatomical volumes for each subject
nconditions = nsessions; % treats each session as a different condition (comment the following three lines and lines 84-86 below if you do not wish to analyze between-session differences)

if nconditions == 1
    batch.Setup.conditions.names = {'rest'};
    for ncond = 1
        for nsub = 1:NSUBJECTS
            for nses = 1:nsessions
                batch.Setup.conditions.onsets{ncond}{nsub}{nses} = 0;
                batch.Setup.conditions.durations{ncond}{nsub}{nses} = inf;
            end
        end
    end % rest condition (all sessions)
else
    batch.Setup.conditions.names = [{'rest'}, arrayfun(@(n)sprintf('Session%d', n), 1:nconditions, 'uni', 0)];
    
    for ncond = 1
        for nsub = 1:NSUBJECTS
            for nses = 1:nsessions
                batch.Setup.conditions.onsets{ncond}{nsub}{nses} = 0;
                batch.Setup.conditions.durations{ncond}{nsub}{nses} = inf;
            end
        end
    end % rest condition (all sessions)
    
    for ncond = 1:nconditions
        for nsub = 1:NSUBJECTS
            for nses = 1:nsessions
                batch.Setup.conditions.onsets{1 + ncond}{nsub}{nses} = [];
                batch.Setup.conditions.durations{1 + ncond}{nsub}{nses} = [];
            end
        end
    end
    
    for ncond = 1:nconditions
        for nsub = 1:NSUBJECTS
            for nses = ncond
                batch.Setup.conditions.onsets{1 + ncond}{nsub}{nses} = 0;
                batch.Setup.conditions.durations{1 + ncond}{nsub}{nses} = inf;
            end
        end
    end % session-specific conditions
end

% Setup ROIs (all anatomical regions, without networks -- Harvard-Oxford atlas)
batch.Setup.rois.files = {fullfile(conn_folder, 'rois', 'atlas.nii')};
load(fullfile(data_dir, 'fcFMRI_analysis', 'ROIs')); % load the ROIs I saved from conn_folder = '/home/uni10/nmri/tools/conn/conn_19b';
batch.Setup.rois.names = ROIs.Regions; % just call the struct field, since I saved it as CELL before
batch.Setup.rois.networks = 0;  % Ensure to turn off network-level analysis

batch.Setup.preprocessing.steps='default_mni'; % RUNS DEFAULT MNI PROCESSING, AVAILABLE IN CONN
% batch.Setup.preprocessing.steps = setdiff(batch.Setup.preprocessing.steps, {'functional_realign&unwarp', 'functional_slice-timing'}); % Remove Realignment & Slice-timing correction steps (SINCE DATA ALREADY PREPROCESSED WITH TOPUP+MCFLIRT WAS SELECTED)
batch.Setup.preprocessing.sliceorder='interleaved (Siemens)';
batch.Setup.done=1;
batch.Setup.overwrite='Yes';

% Run one step at a time:
% conn_batch(batch); % runs Preprocessing and Setup steps only
% clear batch;
% batch.filename=fullfile(cwd,'rs_regions.mat');            % Existing conn_*.mat experiment name

%% DENOISING step
% CONN Denoising                        % Default options (uses White Matter+CSF+realignment+scrubbing+conditions as confound regressors); see conn_batch for additional options
batch.Denoising.filter=[0.008, 0.1];    % frequency filter (band-pass values, in Hz)
                                        % best for resting-state: low-frequency focused filter (0.008-0.1 Hz)
                                        % best for task-based: a high-pass filter (0.008-0.2) keeps the higher-frequency information, more related to BOLD during the events/blocks (more frequent fluctuations)
batch.Denoising.done=1;
batch.Denoising.overwrite='Yes';

% Run one step at a time:
% conn_batch(batch); % runs Denoising step only
% clear batch;
% batch.filename=fullfile(cwd,'rs_regions.mat');            % Existing conn_*.mat experiment name

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
conn('load',fullfile(cwd,'rs_regions.mat'));
conn gui_results
