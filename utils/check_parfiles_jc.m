% DynaCat + Staticat: Anat ROI RSMs, individual stimuli
% 
%   -Computes multivoxel and plots RSMs using indiviudal stimuli as conditions 
%       for each subject and then across subjects within anatomical ROIs.
%   -Choose an ampType. (ampType = 'zscore' or ampType = 'subtractedBetas')
% 
% BR
% 14 May 2024

%% set parameters
subjects = {'AX', 'BR', 'CT', 'DO', 'HC', 'KP', 'RU', 'VL'};
spd = {};
fd = {};
fpd = {};
labels = {};

% Set parameters
numconds = 357;
eventsPerBlock=4; % length of trial is 4s and the TR is 1s

ExpDir='/share/kalanit/biac2/kgs/projects/DynaCat/';
addpath(genpath('/share/kalanit/biac2/kgs/projects/DynaCat/code/scripts/rsms/'))

all_seen_conds = {};

%% compute and plot mvp and rsm per subject and across subjects
for s=1:length(subjects)
    % ===== get subject session info =====% DynaCat + Staticat: Anat ROI RSMs, individual stimuli
    subject = subjects{s};
    [session_path, list_runs, datatype, scan] = dynacat_staticat_sessions(subject);
    cd(session_path)

    % ===== add individual stim parfiles =====
    parfiles = dir(fullfile(session_path, 'parfiles', 'individual_stim_parfiles_set1/'));
    newscans_parfiles = {parfiles(3:end).folder}';
    newscans_parfiles = fullfile(newscans_parfiles, {parfiles(3:end).name}');

    % ===== find condition info from parfiles =====,
    parfile_data = {};
    seen_conds = [];
    for i=1:length(newscans_parfiles)
        file_path = newscans_parfiles{i};
        fid = fopen(file_path, 'r');
        data = textscan(fid, '%s %f %s %f %f %f', 'Delimiter', '\t');
        fclose(fid);
        for ii=1:length(data{4})
            color_1 = data{4}(ii);
            color_2 = data{5}(ii);
            color_3 = data{6}(ii);
            color_code = [color_1 color_2 color_3];
            onset = str2double(data{1}(ii));
            cond_num = data{2}(ii);
            cond = data{3}(ii);
            line_data = [i onset cond_num cond color_code];
            parfile_data = [parfile_data; line_data];
            seen_conds = union(seen_conds, cond_num);
        end
    end
    sorted_parfile_data = sortrows(parfile_data, 3); % sort parfile data based on condition number
    spd{s} = sorted_parfile_data;
    all_seen_conds = [all_seen_conds; seen_conds];

    % remove baseline and task
    filtered_data = sorted_parfile_data(~strcmp(sorted_parfile_data(:, 4), 'baseline') & ~strcmp(sorted_parfile_data(:, 4), 'task-repeat'), :);
    fd{s} = filtered_data;
    
    % remove some repeating condition labels that snuck their way in
    unique_items = {};
    rows_to_keep = true(size(filtered_data, 1), 1);
    for i = 1:size(filtered_data, 1)
        item = filtered_data{i, 4};
        if ismember(item, unique_items)
            rows_to_keep(i) = false;
        else
            unique_items = [unique_items; item];
        end
    end
    filtered_parfile_data = filtered_data(rows_to_keep, :);
    fpd{s} = filtered_parfile_data;

    labels{s} = filtered_parfile_data(:, 4)';
end

br_miss = [];
for i = 1:length(br)-1
    if br(i) < br(i + 1) - 1
        for j = br(i)+1:br(i+1)-1
            br_miss = [br_miss; j];
        end
    end
end