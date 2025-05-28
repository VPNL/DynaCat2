%% DynaCat + StatiCat: Distinctiveness
% Compute within-format and across-format distinctiveness using
% the provided util
%
% JC, February 2025


%% Setup
ExpDir = '/share/kalanit/biac2/kgs/projects/DynaCat/';
ImageDir = fullfile(ExpDir, 'results', 'anat_rois', 'scoring');
DataDir = fullfile(ExpDir, 'results', 'distinctiveness');

addpath(genpath(fullfile(ExpDir, 'code', 'scripts', 'rsms')));
subjects = {'AX', 'BR', 'CT', 'DO', 'HC', 'IK', 'KP', 'RU', 'RY', 'VL'};
formats = {'category', 'animacy', 'motion'}; % 'identity'

ROI_types = {'IF_DO'};
ROIs = cell(1, length(ROI_types)*2);
for r=1:length(ROI_types)
    ROIs{r*2 - 1} = strcat('lh_', ROI_types{r});
    ROIs{r*2} = strcat('rh_', ROI_types{r});
end

%% GLM parameters
eventsPerBlock = 4; % length of trial is 4s and the TR is 1s
ampType = 'zscore';
plot_type = 'sym_means';
ve = 20;

num_all_conds = 384;
parfile_folders = 'individual_stim_parfiles_set1';

category_labels = {'words-d', 'dogs-d', 'bodies-d', 'hands-d', 'faces-d', 'cars-d', 'balls-d', 'scenes-d'...
    'words-s', 'dogs-s', 'bodies-s', 'hands-s', 'faces-s', 'cars-s', 'balls-s', 'scenes-s'};

%% MDS info
stimuli_icon_folder = '/share/kalanit/biac2/kgs/projects/DynaCat/stimuli/stimuli_icons/stimuli_icons_lower_res/dynacat_staticat_colored_icons';
stimuli_icon_label = 'exp_colors';
stimuli_icon_size = 0.015;
categeory_color_point_size = 17;

% load stimuli condition info
all_stimuli_cond_parfile_path = fullfile(ExpDir, 'code', 'experiment_code', 'DynaCat', 'dynacat_stimuli', 'phase_scrambled_set1', 'all_stimuli_conditions.par');
all_cond_parfile = {};
fid = fopen(all_stimuli_cond_parfile_path, 'r');
data = textscan(fid, '%f %s %s %s %f %f %f', 'Delimiter', '\t');length(colors)
fclose(fid);
for i=1:length(data{1})
    cond_num = data{1}(i);
    cond = data{2}(i);
    experiment = data{3}(i);
    category = data{4}(i);
    color_1 = data{5}(i);
    color_2 = data{6}(i);
    color_3 = data{7}(i);
    color_code = [color_1 color_2 color_3];
    line_data = [cond_num cond experiment category color_code];
    all_cond_parfile = [all_cond_parfile; line_data];
end
all_cond_labels = all_cond_parfile(:, 2)';
all_cond_label_colors = all_cond_parfile(:, 5);

% make an array of point shapes for the mds
all_cond_point_shapes = cell(size(all_cond_labels));
for i = 1:numel(all_cond_labels)
    if contains(all_cond_labels{i}, 'dynacat')
        all_cond_point_shapes{i} = 'o';
    elseif contains(all_cond_labels{i}, 'staticat')
        all_cond_point_shapes{i} = '^';
    end
end

%% Compute subject-level MVPs + RSMs
for r = 1:length(ROIs)
    ROIname = ROIs{r};
    fprintf('Generating %s RSM at %.2f pct VE\n\n', ROIname, ve);
    %fprintf('Generating %s RSM\n\n', ROIname);

    % initialize across-subject rsm + distance matrix
    across_subject_rsm = nan(num_all_conds, num_all_conds, length(subjects));

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

        % get betas for each run and voxel
        amps = mv_amps(mv);

        % removes task condition from amps
        amps = amps(:, 1:numconds);

        % values to plot
        z_values = amps;

        if strcmp(ampType, 'zscore')
            % subtract mean beta from each voxel
            mean_amps=mean(amps,2); % calculate the mean for each voxel
            meanmat=mean_amps*ones(1,numconds);
            subbetas=amps-meanmat;
            % normalize by the residual variance of the GLM in each voxel
            % effectively after this step you have z-scoed your pattern
            residualV = mv.glm.residual; % make this into functions
            dof = mv.glm.dof;
            residualvar = sum(residualV.^2)/dof; % calculate mean sqaure error
            resd = sqrt(residualvar);
            resd = resd';
            resd_mat = repmat(resd, [1 numconds]); %turn into matrix
            z_values = subbetas./resd_mat; % divide each beta by standard deviation of residual error
        end

        % titrate
        voxels_ve = find(mv.glm.varianceExplained > ve/100);
        z_values = z_values(voxels_ve',:);

        % compute crosscorrelation matrices
        rsm = zeros(numconds);
        all_cond_rsm = NaN(num_all_conds, num_all_conds);
        for i=1:numconds % calculate correlation between each condition
            cond_num_x = cell2mat(filtered_parfile_data(i, 3));
            for j=1:numconds
                c = corrcoef(z_values(:,i),z_values(:,j));
                rsm(i, j) = c(1, 2); % subject rsm

                cond_num_y = cell2mat(filtered_parfile_data(j, 3));
                all_cond_rsm(cond_num_x, cond_num_y) = rsm(i, j); % all cond subject rsm
            end
        end

        % save subject rsm data
        subjectDataDir = fullfile(DataDir, subject);
        if ~exist(subjectDataDir, 'dir')
            mkdir(subjectDataDir);
        end
        rsm_filepath = fullfile(subjectDataDir, [ROIname '_' num2str(ve) '_ve_rsm.mat']);
        %rsm_filepath = fullfile(subjectDataDir, [ROIname '_rsm.mat']);
        save(rsm_filepath, 'rsm');
        % save subject all cond rsm data
        rsm_all_cond_filepath = fullfile(subjectDataDir, [ROIname '_' num2str(ve) '_ve_rsm_all_cond.mat']);
        %rsm_all_cond_filepath = fullfile(subjectDataDir, [ROIname '_rsm_all_cond.mat']);
        save(rsm_all_cond_filepath, 'all_cond_rsm');

        across_subject_rsm(:, :, s) = all_cond_rsm; % add this subject to across subject rsm

        % clear, close everything to save memory
        close all
        clear roi_rsm
        clear rsm
        clear all_cond_rsm
        clear mv
        clear z_values
        clear amps
        clear allvals
        clear cmap
        clear mean_amps
        clear meanmat
        clear mvs
        clear resd
        clear resd_mat
        clear residualV
        clear residualvar
        clear subbetas
    end

    % get mean across subjects rsms per roi
    mean_rsm = nanmean(across_subject_rsm, 3);

    % save mean rsm data
    mean_rsm_filepath = fullfile(DataDir, [ROIname '_' num2str(ve) '_ve_rsm.mat']);
    %mean_rsm_filepath = fullfile(DataDir, [ROIname '_rsm.mat']);
    save(mean_rsm_filepath, 'mean_rsm');

    close all 
    clear mean_rsm
end


%% Compute distinctiveness
% store each distinctiveness score in a 3D matrix, indexed score type x ROI x
% subject
fprintf('Computing distinctiveness...\n');

n_scoretypes = 7;
scores = nan(n_scoretypes, length(ROIs), length(subjects));
for s = 1:length(subjects)
    subject = subjects{s};
    for r = 1:length(ROIs)
        ROI = ROIs{r};
        filename = [ROI '_rsm_all_cond.mat'];
        load(fullfile(DataDir, subject, filename));

        scores(1:3, r, s) = distinctiveness_category(all_cond_rsm);
        scores(4:6, r, s) = distinctiveness_animacy(all_cond_rsm);
        scores(7, r, s) = distinctiveness_format(all_cond_rsm);
    end
end

cd(DataDir);
save('scores.mat', 'scores', '-append');
fprintf('Scores saved and calculated!\n\n');

%% Plotting distinctiveness
fprintf('Plotting distinctiveness...\n');
avg_scores = mean(scores, 3);
std_scores = std(scores, 0, 3);
score_types = {'cat (dynamic)', 'cat (static)', 'cat (across)', ...
    'anim (dynamic)', 'anim (static)', 'anim (across)', 'format'};
FigDir = fullfile(DataDir, 'figures');
cd(FigDir);

for r = 1:length(ROIs)
    ROI = ROIs{r};
    plotDistinctiveness(ROI, score_types, avg_scores(:,r), std_scores(:,r));
end
fprintf('All done plotting!\n\n');


%% Distinctiveness per category
% store each distinctiveness in a 3D matrix, indexed category (1-8) x ROI x
% subject
fprintf('Calculating distinctiveness per category...\n');

n_categories = 8;
cat_scores_across = nan(n_categories, length(ROIs), length(subjects));
cat_scores_d = nan(n_categories, length(ROIs), length(subjects));
cat_scores_s = nan(n_categories, length(ROIs), length(subjects));

for s = 1:length(subjects)
    subject = subjects{s};
    for r = 1:length(ROIs)
        ROI = ROIs{r};
        filename = [ROI '_rsm_all_cond.mat'];
        load(fullfile(DataDir, subject, filename));

        all_scores = distinctiveness_per_category(all_cond_rsm);
        cat_scores_across(1:8, r, s) = all_scores(1,:,:);
        cat_scores_d(1:8, r, s) = all_scores(2,:,:);
        cat_scores_s(1:8, r, s) = all_scores(3,:,:);
    end
end

cd(DataDir);
save('scores.mat', 'cat_scores_across', 'cat_scores_d', 'cat_scores_s', '-append');
fprintf('Done calculating distinctiveness per category!\n\n');

%% Plotting distinctiveness per category
fprintf('Plotting distinctiveness per category...\n');
category_names = {'words', 'dogs', 'bodies', 'hands', 'faces', 'cars', ...
    'balls', 'scenes'};
FigDir = fullfile(DataDir, 'figures', 'per_cat');
cd(FigDir);

for r = 1:length(ROIs)
    ROI = ROIs{r};
    scores_across = cat_scores_across(:,r,:);
    scores_d = cat_scores_d(:,r,:);
    scores_s = cat_scores_s(:,r,:);
    means = [mean(scores_d, 3) mean(scores_s, 3) mean(scores_across, 3)];
    stds = [std(scores_d, 0, 3) std(scores_s, 0, 3) std(scores_across, 0, 3)];
    plotDistinctivenessPerCat(ROI, category_names, means, stds);
end
fprintf('All done plotting!\n\n');

%% Functions
function [] = plotDistinctiveness(ROI, xlabs, data, errors)
hBar = bar(data, 'facecolor', 'flat');
errhigh = [];
errlow = [];
for i = 1:length(data)
    errhigh = [errhigh (data(i) - errors(i))];
    errlow = [errlow (data(i) - errors(i))];
end
ylim([-0.01 0.15]);

hold on
er = errorbar(1:length(data), data, errlow, errhigh);
er.Color = [0 0 0];
er.LineStyle = 'none';

clr = [0 0.447 0.741; 0 0.447 0.741; 0 0.447 0.741; 0.85 0.325 0.098; ...
    0.85 0.325 0.098; 0.85 0.325 0.098; 0.929 0.694 0.125];
hBar.CData = clr;

title([ROI ' Distinctiveness Scores']);
ax = hBar.Parent;
set(ax, 'xtick', 1:length(data));
set(gca, 'xticklabel', xlabs);
set(gcf, 'position', [100, 100, 1700, 1000]);

filename = [ROI '_distinctiveness_plot.png'];
saveas(gcf, filename);
hold off
close
end

function [] = plotDistinctivenessPerCat(ROI, xlabs, data, errors)
hBar = bar(data, 'facecolor', 'flat');
hold on

ylim([-0.03 0.24]);

ngroups = size(data, 1);
nbars = size(data, 2);
groupwidth = min(0.8, nbars / (nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i - 1) * groupwidth / (2 * nbars);
    er = errorbar(x, data(:, i), errors(:, i), '.');
    er.Color = [0 0 0];
end

title([ROI ' Distinctiveness Per Category']);
ax = hBar.Parent;
set(ax, 'xtick', 1:length(data));
set(gca, 'xticklabel', xlabs);
set(gcf, 'position', [100, 100, 1700, 1000]);

filename = [ROI '_distinctiveness_per_cat_plot.png'];
saveas(gcf, filename);
hold off
close
end

% returns elements above the diagonal of a square matrix as an array
function [flattened] = aboveDiagonalFlat(M)
flattened = [];
for row=1:length(M)
    flattened = [flattened M(row, row+1:end)];
end
end


% returns elements of a matrix as an array
function [flattened] = flat(M)
flattened = reshape(M.',1,[]);
end

% dynamic, static, both
function [scores] = distinctiveness_per_category(rsm)
cat_scores = zeros(1, 8);
cat_scores_d = zeros(1, 8);
cat_scores_s = zeros(1, 8);

for condition=1:8
    start_dyn = ((condition - 1) * 24) + 1;
    start_sta = start_dyn + 192;

    % without cluster
    dyn_to_other = rsm(start_dyn:start_dyn + 23, [1:start_dyn - 1, ...
        start_dyn + 24:192]);
    sta_to_other = rsm(start_sta:start_sta + 23, [193:start_sta - 1, ...
        start_sta + 24:end]);
    dyn_to_sta_other = rsm(start_dyn:start_dyn + 23, [193:start_sta - 1, ...
        start_sta + 24:end]);
    sta_to_dyn_other = rsm(start_sta:start_sta + 23, [1:start_dyn - 1, ...
        start_dyn + 24:192]);
    without_cluster = [flat(dyn_to_other) flat(sta_to_other) flat(dyn_to_sta_other) ...
        flat(sta_to_dyn_other)];
    a = nanmean(without_cluster);
    a_d = nanmean(flat(dyn_to_other));
    a_s = nanmean(flat(sta_to_other));

    % within cluster
    dyn_cluster = aboveDiagonalFlat(rsm(start_dyn:start_dyn + 23, start_dyn:start_dyn + 23));
    sta_cluster = aboveDiagonalFlat(rsm(start_sta:start_sta + 23, start_sta:start_sta + 23));
    across_format_within_cat = rsm(start_dyn:start_dyn + 23, start_sta:start_sta + 23);
    within_cluster = [dyn_cluster sta_cluster flat(across_format_within_cat)];
    b = nanmean(within_cluster);
    b_d = nanmean(dyn_cluster);
    b_s = nanmean(sta_cluster);

    cat_scores(condition) = (b - a);
    cat_scores_d(condition) = (b_d - a_d);
    cat_scores_s(condition) = (b_s - a_s);
end
scores = [cat_scores; cat_scores_d; cat_scores_s];
%score_b = mean(cat_scores);
%score_d = mean(cat_scores_d);
%score_s = mean(cat_scores_s);

%scores = [score_d score_s score_b];
end

% dynamic, static, both
function [scores] = distinctiveness_animacy(rsm)
% across format
ani_to_ani_dyn = aboveDiagonalFlat(rsm(25:120, 25:120));
ani_to_ani_sta = aboveDiagonalFlat(rsm(217:312, 217:312));
ani_to_ani_across = rsm(25:120, 217:312);
within_cluster_ani = [ani_to_ani_dyn ani_to_ani_sta ...
    flat(ani_to_ani_across)];
b_ani = nanmean(within_cluster_ani);

ina_words = aboveDiagonalFlat(rsm(1:24, 1:24));
ina_mixed = aboveDiagonalFlat(rsm(121:216, 121:216));
ina_sta = aboveDiagonalFlat(rsm(313:end, 313:end));
within_cluster_ina = [ina_words ina_mixed ina_sta ...
    flat(rsm(1:24, [121:216, 313:end])) flat(rsm(121:216, 313:end))];
b_ina = nanmean(within_cluster_ina);

without_cluster = [flat(rsm(1:24, [25:120, 217:312])) ...
    flat(rsm(25:120, [121:216, 313:end])) ...
    flat(rsm(121:216, 217:312)) ...
    flat(rsm(217:312, 313:end))];
a = nanmean(without_cluster);

distinctiveness_ani = (b_ani - a);
distinctiveness_ina = (b_ina - a);
score_b = mean([distinctiveness_ani distinctiveness_ina]);

% dynamic
b_ani_d = nanmean(ani_to_ani_dyn);
b_ina_d = nanmean([ina_words, flat(rsm(1:24, 121:192)), ...
    aboveDiagonalFlat(rsm(121:192, 121:192))]);
a_d = nanmean([flat(rsm(1:24, 25:120)), flat(rsm(25:120, 121:192))]);
score_d = mean([(b_ani_d - a_d) (b_ina_d - a_d)]);

% static
b_ani_s = nanmean(ani_to_ani_sta);
b_ina_s = nanmean([aboveDiagonalFlat(rsm(193:216, 193:216)), flat(rsm(193:216, 313:end)), ...
    aboveDiagonalFlat(rsm(313:end, 313:end))]);
a_s = nanmean([flat(rsm(193:216,217:312)), flat(rsm(217:312, 313:end))]);
score_s = mean([(b_ani_s - a_s) (b_ina_s - a_s)]);

scores = [score_d score_s score_b];
end

function [score] = distinctiveness_format(rsm)
dyn_to_dyn = nanmean(aboveDiagonalFlat(rsm(1:192, 1:192)));
sta_to_sta = nanmean(aboveDiagonalFlat(rsm(193:end, 193:end)));
dyn_to_sta = nanmean(flat(rsm(1:192, 193:end)));

distinctiveness_dyn = dyn_to_dyn - dyn_to_sta;
distinctiveness_sta = sta_to_sta - dyn_to_sta;

score = mean([distinctiveness_dyn distinctiveness_sta]);
end