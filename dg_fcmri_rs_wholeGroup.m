% dg_fcmri_rs_draft
clear all; close all; clc;
cd('PROJECT/DIRECTORY');
% run('dg_open_Conn_19b.m'); % LOADS SPM12 AND CONN
addpath(genpath('updateThisDirectory/spm12_7487'));
addpath(genpath('updateThisDirectory/conn_19b'));

%% FIND AND ADD functional/structural files
% note: this will look for all data in these folders, irrespective of the specific download subsets entered as command-line arguments
data_dir = 'PROJECT/DIRECTORY';
cwd = 'RESULTS/DIRECTORY/fcFMRI_analysis/rs';

patients_list = {'pseudonyms'};
controls_list = {'pseudonyms'};

% patient's functional files (UNPROCESSED)
pf_files = [];
for i = 1:length(patients_list)
    pf_files{i} = fullfile(data_dir, 'ROLANDIC', patients_list{i}, 'func', '--dataset--', [patients_list{i}, '.EPI_fMRI.rs.PA.nii.gz']);
end
pf_files = pf_files';

% control's functional files (UNPROCESSED)
cf_files = [];
for i = 1:length(controls_list)
    cf_files{i} = fullfile(data_dir, 'CONTROLS', controls_list{i}, 'func', '--dataset--', [controls_list{i}, '.EPI_fMRI.rs.PA.nii.gz']);
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
batch.filename=fullfile(cwd,'rs.mat');            % New conn_*.mat experiment name

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

batch.Setup.preprocessing.steps='default_mni'; % RUNS DEFAULT MNI PROCESSING, AVAILABLE IN CONN
batch.Setup.preprocessing.sliceorder='interleaved (Siemens)';
batch.Setup.done=1;
batch.Setup.overwrite='Yes';

% Run one step at a time:
% conn_batch(batch); % runs Preprocessing and Setup steps only
% clear batch;
% batch.filename=fullfile(cwd,'rs.mat');            % Existing conn_*.mat experiment name

%% DENOISING step
% CONN Denoising                                    % Default options (uses White Matter+CSF+realignment+scrubbing+conditions as confound regressors); see conn_batch for additional options
batch.Denoising.filter=[0.01, 0.1];                 % frequency filter (band-pass values, in Hz)
batch.Denoising.done=1;
batch.Denoising.overwrite='Yes';

% Run one step at a time:
% conn_batch(batch); % runs Denoising step only
% clear batch;
% batch.filename=fullfile(cwd,'rs.mat');            % Existing conn_*.mat experiment name

%% FIRST-LEVEL ANALYSIS step
% CONN Analysis                                     % Default options (uses all ROIs in conn/rois/ as connectivity sources); see conn_batch for additional options
batch.Analysis.done=1;
batch.Analysis.overwrite='Yes';

% Run ALL analyses
conn_batch(batch);

% CONN Display
% launches conn gui to explore results AND RUN 2nd LEVEL ANALYSIS WITH GUI
conn
conn('load',fullfile(cwd,'rs.mat'));
conn gui_results