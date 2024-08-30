function [subject] = dg_child_taskdata_trim(subject)
% This script extracts all task and scanner timepoints, scrubbing data, and movement parameters to stablish regressors in the fMRI model
%
% KEY -> The current subject structure must be loaded to specify folders and extract all the variables content!
%
% A TRIMMING threshold must be set:
trimming_threshold = 0.1; % In case of datasets longer than the task duration, TRIM if the difference is > 10% of the total dataset length

if contains(subject.exam_id, 'language')
    clearvars -except subject trimming_threshold
    condition_folder = 'wl';
    fprintf('Extracting LANGUAGE task and scan data for subject %s \n', subject.id);
    
    %% Extract task timepoints
    filename = [subject.analysis_dir '/' subject.id '/' condition_folder '/' 'wl.txt' ];
    
    delimiter = '\t';
    startRow = 2;
    formatSpec = '%s%s%f%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    wl = table(dataArray{1:end-1}, 'VariableNames', {'imagesDir','image_shownname','Accuracy','Reaction_time','image_shown','Accumulated_time'});
    clearvars filename delimiter startRow formatSpec fileID dataArray ans;
    
    % Animals and Fractals tables by indexing (INDEXING ACCORDING TO COUNTERBALLANCING)
    text_images = 'mixed_images';
    text_fractals = 'mixed_fractals';
    
    for i = 1:size(wl,1)
        cond = string(wl{i,1});
        trials_images(i) = strcmp(cond, text_images);
        trials_fractals(i)  = strcmp(cond, text_fractals);
    end
    
    Animals = wl(trials_images,:);
    Fractals = wl(trials_fractals,:);
    
    % GET ONSETS (TIMEPOINTS) (FROM ACCUMULATED TIME)
    subject.WLtask.Onsets_animals = table2array(Animals(:,6));
    subject.WLtask.Onsets_fractals = table2array(Fractals(:,6));
    
    % GET ACCURACY AND REACTION TIMES
    subject.WLtask.Accuracy_animals = table2array(Animals(:,3));
    subject.WLtask.Accuracy_fractals = table2array(Fractals(:,3));
    
    subject.WLtask.RT_animals = table2array(Animals(:,4));
    subject.WLtask.RT_fractals = table2array(Fractals(:,4));
    
    %% Extract movement parameters (movement paramaters from .par file, without ART scrubbing!)
    filename = [subject.analysis_dir '/' subject.id '/func/' subject.id '.EPI_fMRI.language.PA.mcf.par'];
    
    delimiter = ' ';
    formatSpec = '%f%f%f%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string',  'ReturnOnError', false);
    fclose(fileID);
    R = [dataArray{1:end-1}];
    clearvars filename delimiter formatSpec fileID dataArray ans;
    
    original_scans = length(R);                         % total number of volumes (from the total number of data in the .par file)
    subject.R_movement_params = R;                      % saves the 6-movement parameters in subject structure
    
    %% Comparing scanner-task timepoints (for eventual trimming of the dataset length - TAKING TASK TIMEPOINTS, not scanner triggers)
    TR = 1.400;                                                             % TR
    task_timepoints = table2array(wl(:,6));                                 % takes the Accumulated times of the Stimuli presentation as task timepoints (both from Animals and Fractals)
    last_timepoint = task_timepoints(end,1);
    last_reaction_time = table2array(wl(end,4));                            % takes the last Reaction Time (RT)
    last_task_timepoint = last_timepoint + last_reaction_time;              % takes the last Accumulated time PLUS the last of the Reaction Times (RT) as the last task timepoint
    
    total_time = original_scans * TR;                                       % calculates the total time of the scan (in seconds)
    scanner_timepoints = linspace(0, total_time, original_scans);           % creates a theoretical timeline of the scanner pulses (taking the total number of volumes and TR)
    
    if last_reaction_time < TR                                              % to add the last RT (if last RT > TR) or add one just one TR (if last RT < TR)
        disp('add the last RT (if last RT > TR).')
        last_task_point = last_task_timepoint + TR;
    elseif last_reaction_time > TR
        disp('add one just one TR (if last RT < TR).')
        last_task_point = last_task_timepoint + last_reaction_time;
    end
    
    % DETERMINE TRIMMING
    % Determine which of the volumes is the closest to the last task timepoint (in this case, last reaction time RT)
    [time_dif, last_task_scan] = min(abs(last_task_point - scanner_timepoints));
    %last_task_scan = interp1(scanner_timepoints, 1:length(scanner_timepoints), last_task_point, 'nearest');
    
    if scanner_timepoints(last_task_scan) <= last_task_point
        disp('Last scanner timepoint closest to last task point is smaller. Using an additional scan to the end.')
        last_task_scan = last_task_scan + 1 ;
    else
        disp('Last scanner timepoint closest to last task point is greater. Found last scan point is correct.')
        last_task_scan = last_task_scan;
    end
    
    difference_volumes = original_scans - last_task_scan;
    comparison_threshold = trimming_threshold * original_scans;
    if original_scans > last_task_scan && difference_volumes > comparison_threshold % 2 criteria: trim if original volumes > task scans AND difference > 10% of original volumes
        fprintf('Scan length trimmed!, difference between total of scans and last task scan > than trimming threshold \n');
        subject.scan_trimmed_length = (1:last_task_scan)';
    elseif (1 < difference_volumes) && (difference_volumes < comparison_threshold)  % Don't trim if difference is between 1 and 10% of original volumes
        fprintf('Scan length not trimmed, difference between total of scans and last task scan < than trimming threshold \n');
        subject.scan_trimmed_length = [];
        subject.scan_length = original_scans;
    elseif original_scans == last_task_scan                                         % Don't trim if there's no difference between original volumes and last task scan
        fprintf('Scan length not trimmed, original total of scans = to last task scan \n');
        subject.scan_trimmed_length = [];
        subject.scan_length = original_scans;
    else
        fprintf('Scan length not trimmed, difference between total of scans and last task scan < than trimming threshold \n');
        subject.scan_trimmed_length = [];
        subject.scan_length = original_scans;
    end
    
    % DECISIONS TO DETERMINE REGRESSORS!
    if ~isempty(subject.scan_trimmed_length)    % IF TRIMMING NEEDED
        subject.R_movement_params = subject.R_movement_params(1:length(subject.scan_trimmed_length),:);
        R = subject.R_movement_params;                                                               % Re-writes R!
        movparams_file = fullfile(subject.analysis_dir, subject.id, '/func/', 'mov_params_wl.mat');
        save(movparams_file', 'R');                                                                  % Saves the TRIMMED 6-movement parameters as the only regressor
    elseif isempty(subject.scan_trimmed_length) % IF NO TRIMMING NEEDED
        R = subject.R_movement_params;                                                               % Re-writes R! Although in this case it is still the same R from the .par file
        movparams_file = fullfile(subject.analysis_dir, subject.id, '/func/', 'mov_params_wl.mat');
        save(movparams_file', 'R');                                                                  % Saves the UNTRIMMED 6-movement parameters as the only regressor
    end
    
    %% Update the subject struct!
    %     params_file = fullfile(subject.analysis_dir, ['analysis_params.m']);
    %     run(params_file);
    dv_update_subjectinfo(subject, subject.analysis_dir);
    
elseif contains(subject.exam_id, 'WM')
    clearvars -except subject trimming_threshold
    condition_folder = 'wm';
    fprintf('Extracting WORKING MEMORY task and scan data for subject %s \n', subject.id);
    
    %% Extract task timepoints
    filename = [subject.analysis_dir '/' subject.id '/' condition_folder '/' 'wm.txt' ];
    
    delimiter = '\t';
    startRow = 2;
    formatSpec = '%f%f%C%C%f%f%s%f%f%f%f%f%f%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    wm = table(dataArray{1:end-1}, 'VariableNames', {'subid','subage','gender','group','blockNumber','trialNumber','presentation','is_correct','keypressed','accuracy','ReactionTime','AccumulatedTimeResponse','AccumulatedTimeEncoding1','AccumulatedTimeEncoding2','AccumulatedTimeEncoding3','AccumulatedTimeRetrieval'});
    clearvars filename delimiter startRow formatSpec fileID dataArray ans;
    
    % GET ONSETS (TIMEPOINTS) (ACCUMULATED TIMES (ENC + RET))
    enc1 = table2array(wm(:,13)); enc2 = table2array(wm(:,14)); enc3 = table2array(wm(:,15));
    enc_not_ordered = [enc1; enc2; enc3];
    
    subject.WMtask.Onsets_encoding = sort(enc_not_ordered, 'ascend');
    subject.WMtask.Onsets_retrieval = table2array(wm(:,16));
    
    % GET ACCURACY AND REACTION TIMES
    subject.WMtask.Accuracy_retrieval = table2array(wm(:,10));
    subject.WMtask.RT_retrieval = table2array(wm(:,11));
    
    %% Extract movement parameters (movement paramaters from .par file, without ART scrubbing!)
    filename = [subject.analysis_dir '/' subject.id '/func/' subject.id '.EPI_fMRI.WM.PA.mcf.par'];
    
    delimiter = ' ';
    formatSpec = '%f%f%f%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string',  'ReturnOnError', false);
    fclose(fileID);
    R = [dataArray{1:end-1}];
    clearvars filename delimiter formatSpec fileID dataArray ans;
    
    original_scans = length(R);                         % total number of volumes (from the total number of data in the .par file)
    subject.R_movement_params = R;                      % saves the 6-movement parameters in subject structure
    
    %% Comparing scanner-task timepoints (for eventual trimming of the dataset length - TAKING TASK TIMEPOINTS, not scanner triggers)
    TR = 1.400;                                                             % TR
    task_timepoints = table2array(wm(:,12));                                % takes the Accumulated times at RESPONSES (responses of retrieval)
    last_timepoint = task_timepoints(end,1);
    last_reaction_time = table2array(wm(end,11));                           % takes the last Reaction Time (RT)
    last_task_timepoint = last_timepoint + last_reaction_time;              % takes the last Accumulated time PLUS the last of the Reaction Times (RT) as the last task timepoint
    
    total_time = original_scans * TR;                                       % calculates the total time of the scan (in seconds)
    scanner_timepoints = linspace(0, total_time, original_scans);           % creates a theoretical timeline of the scanner pulses (taking the total number of volumes and TR)
    
    if last_reaction_time < TR                                              % to add the last RT (if last RT > TR) or add one just one TR (if last RT < TR)
        disp('add the last RT (if last RT > TR).')
        last_task_point = last_task_timepoint + TR;
    elseif last_reaction_time > TR
        disp('add one just one TR (if last RT < TR).')
        last_task_point = last_task_timepoint + last_reaction_time;
    end
    
    % DETERMINE TRIMMING
    % Determine which of the volumes is the closest to the last task timepoint (in this case, last reaction time RT)
    [time_dif, last_task_scan] = min(abs(last_task_point - scanner_timepoints));
    %last_task_scan = interp1(scanner_timepoints, 1:length(scanner_timepoints), last_task_point, 'nearest');
    
    if scanner_timepoints(last_task_scan) <= last_task_point
        disp('Last scanner timepoint closest to last task point is smaller. Using an additional scan to the end.')
        last_task_scan = last_task_scan + 1 ;
    else
        disp('Last scanner timepoint closest to last task point is greater. Found last scan point is correct.')
        last_task_scan = last_task_scan;
    end
    
    difference_volumes = original_scans - last_task_scan;
    comparison_threshold = trimming_threshold * original_scans;
    if original_scans > last_task_scan && difference_volumes > comparison_threshold % 2 criteria: trim if original volumes > task scans AND difference > 10% of original volumes
        fprintf('Scan length trimmed!, difference between total of scans and last task scan > than trimming threshold \n');
        subject.scan_trimmed_length = (1:last_task_scan)';
    elseif (1 < difference_volumes) && (difference_volumes < comparison_threshold)  % Don't trim if difference is between 1 and 10% of original volumes
        fprintf('Scan length not trimmed, difference between total of scans and last task scan < than trimming threshold \n');
        subject.scan_trimmed_length = [];
        subject.scan_length = original_scans;
    elseif original_scans == last_task_scan                                         % Don't trim if there's no difference between original volumes and last task scan
        fprintf('Scan length not trimmed, original total of scans = to last task scan \n');
        subject.scan_trimmed_length = [];
        subject.scan_length = original_scans;
    else
        fprintf('Scan length not trimmed, difference between total of scans and last task scan < than trimming threshold \n');
        subject.scan_trimmed_length = [];
        subject.scan_length = original_scans;
    end
    
    % DECISIONS TO DETERMINE REGRESSORS!
    if ~isempty(subject.scan_trimmed_length)    % IF TRIMMING NEEDED
        subject.R_movement_params = subject.R_movement_params(1:length(subject.scan_trimmed_length),:);
        R = subject.R_movement_params;                                                               % Re-writes R!
        movparams_file = fullfile(subject.analysis_dir, subject.id, '/func/', 'mov_params_wm.mat');
        save(movparams_file', 'R');                                                                  % Saves the TRIMMED 6-movement parameters as the only regressor
    elseif isempty(subject.scan_trimmed_length) % IF NO TRIMMING NEEDED
        R = subject.R_movement_params;                                                               % Re-writes R! Although in this case it is still the same R from the .par file
        movparams_file = fullfile(subject.analysis_dir, subject.id, '/func/', 'mov_params_wm.mat');
        save(movparams_file', 'R');                                                                  % Saves the UNTRIMMED 6-movement parameters as the only regressor
    end
    
    %% Update new subject struct
    %     params_file = fullfile(subject.analysis_dir, ['analysis_params.m']);
    %     run(params_file);
    dv_update_subjectinfo(subject);
    
elseif contains(subject.exam_id, 'rs')
    clearvars -except subject trimming_threshold
    condition_folder = 'rs';
    fprintf('Extracting RESTING STATE condition and scan data for subject %s \n', subject.id);
    
    % WHATEVER I MAY NEED LATER
    
    %% Update new subject struct
    %     params_file = fullfile(subject.analysis_dir, ['analysis_params.m']);
    %     run(params_file);
    dv_update_subjectinfo(subject);
    
else
    error('No condition identified in raw data and in the experiment design \n')
    
end
end
