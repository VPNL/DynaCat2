%% DynaCat + StatiCat: Mean VE
%
% JC, February 2025


%% Setup
ExpDir = '/share/kalanit/biac2/kgs/projects/DynaCat/';
ImageDir = fullfile(ExpDir, 'results', 'anat_rois', 'scoring');
DataDir = fullfile(ExpDir, 'code', 'justin', 'rsms');

addpath(genpath(fullfile(ExpDir, 'code', 'scripts', 'rsms')));
subjects = {'AX', 'BR', 'CT', 'DO', 'HC', 'IK', 'KP', 'RU', 'RY', 'VL'};
formats = {'category', 'animacy', 'motion'}; % 'identity'

ROI_types = {'IPS' , 'STS', 'LOTC', 'VTC', 'IF_DO'};
ROIs = cell(1, length(ROI_types)*2);
for r=1:length(ROI_types)
    ROIs{r*2 - 1} = strcat('lh_', ROI_types{r});
    ROIs{r*2} = strcat('rh_', ROI_types{r});
end

%% GLM parameters
eventsPerBlock = 4; % length of trial is 4s and the TR is 1s
ampType = 'zscore';
plot_type = 'sym_means';

num_all_conds = 384;
parfile_folders = 'individual_stim_parfiles_set1';

category_labels = {'words-d', 'dogs-d', 'bodies-d', 'hands-d', 'faces-d', 'cars-d', 'balls-d', 'scenes-d'...
    'words-s', 'dogs-s', 'bodies-s', 'hands-s', 'faces-s', 'cars-s', 'balls-s', 'scenes-s'};


%% Compute mean VE
mean_ve_list = nan(length(ROIs), length(subjects));
for r = 1:length(ROIs)
    ROIname = ROIs{r};

    for s = 1:length(subjects)
        % ===== get subject session info =====
        subject = subjects{s};
        [session_path, list_runs, datatype, scan] = dynacat_staticat_sessions(subject);
        cd(session_path)

        % ===== add individual stim parfiles =====
        parfiles = dir(fullfile(session_path, 'parfiles', parfile_folders));
        newscans_parfiles = {parfiles(3:end).folder}';
        newscans_parfiles = fullfile(newscans_parfiles, {parfiles(3:end).name}');

        % ===== find condition info from parfiles =====
        parfile_data = {};
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
            end
        end
        sorted_parfile_data = sortrows(parfile_data, 3); % sort parfile data based on condition number

        % remove baseline and task
        cond_nums = cell2mat(sorted_parfile_data(:, 3));
        rowsToKeep = (cond_nums ~= 0) & (cond_nums ~= 999);
        filtered_data = sorted_parfile_data(rowsToKeep, :);

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
        parfile_cond_nums = cell2mat(filtered_parfile_data(:, 3));
        numconds = length(unique(parfile_cond_nums));

        % ===== plot rsm =====
        % copy ROI over from mrVistaROI folder if nesecary
        ROI_filepath = fullfile(session_path, '3DAnatomy', 'ROIs', [ROIname '.mat']);
        ROI_filepath_in_mrVistaROI_dir = fullfile(session_path, '3DAnatomy', 'mrVistaROIs', [ROIname '.mat']);
        if ~exist(ROI_filepath, 'file')
            copyfile(ROI_filepath_in_mrVistaROI_dir, ROI_filepath);
        end

        % set vANATOMY
        hg = initHiddenGray(datatype, scan, ROIname);

        % add individual stim parfiles
        runs = list_runs;
        hg = er_assignParfilesToScans(hg, runs, newscans_parfiles);

        % initialize multivoxel pattern
        mv = mv_init(hg,ROIname,runs,datatype);
        mvs = {mv};

        % set parameters
        params = mvs{1}.params;
        params.eventAnalysis = 1;
        params.detrend = 1; % high-pass filter (remove the low frequency trend)
        params.detrendFrames = 20;
        params.inhomoCorrection = 1; % divide by mean (transforms from raw scanner units to % signal)
        params.temporalNormalization = 0; % no matching 1st temporal frame
        params.ampType = 'betas';
        params.glmHRF = 3; % SPM hrf
        params.eventsPerBlock = eventsPerBlock;
        params.lowPassFilter = 0; % no temporal low pass filter of the data

        mv.params = params;

        % apply GLM to multivoxel pattern
        mv = mv_applyGlm(mv);

        % store value
        mean_ve_list(r, s) = mean(mv.glm.varianceExplained);
    end
end

%% Plotting mean VE
data = mean(mean_ve_list, 2);
errors = std(mean_ve_list, 0, 2);
FigDir = fullfile(DataDir, 'figures');
cd(FigDir);

plotMeanVE(ROIs, data, errors);

%% Functions
function [] = plotMeanVE(xlabs, data, errors)
hBar = bar(data, 'facecolor', 'flat');
errhigh = [];
errlow = [];
for i = 1:length(data)
    errhigh = [errhigh (data(i) - errors(i))];
    errlow = [errlow (data(i) - errors(i))];
end
ylim([0 1]);

hold on
er = errorbar(1:length(data), data, errlow, errhigh);
er.Color = [0 0 0];
er.LineStyle = 'none';

clr = [0 0.447 0.741; 0.85 0.325 0.098; 0 0.447 0.741; 0.85 0.325 0.098; ...
    0 0.447 0.741; 0.85 0.325 0.098; 0 0.447 0.741; 0.85 0.325 0.098; ...
    0 0.447 0.741; 0.85 0.325 0.098];
hBar.CData = clr;

title('Mean VE Across ROIs');
ax = hBar.Parent;
set(ax, 'xtick', 1:length(data));
set(gca, 'xticklabel', xlabs);
set(gcf, 'position', [100, 100, 1700, 1000]);

filename = 'mean_ve.png';
saveas(gcf, filename);
hold off
close
end
