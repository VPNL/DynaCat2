%% DynaCat + StatiCat: Distinctiveness
% Compute within-format and across-format distinctiveness
% 
% JC, November 2024


%% Setup
ExpDir = '/share/kalanit/biac2/kgs/projects/DynaCat/';
ResultsDir = fullfile(ExpDir, 'results', 'func_rois', 'scoring');
ImageDir = fullfile(ResultsDir, 'Images');
FigDir = fullfile(ResultsDir, 'Figures');
addpath(genpath(fullfile(ExpDir, 'code', 'scripts', 'rsms')));
subjects = {'AX', 'BR', 'CT', 'DO', 'HC', 'IK', 'KP', 'RU', 'RY', 'VL'};

ROI_types = {'pfus_faces' 'mfus_faces' 'ots_bodies' 'los_bodies' 'itg_bodies' 'mtg_bodies' ...
    'iog_faces' 'psts_faces' 'pots_words' 'ios_words' 'ips_scenes' 'mog_scenes' 'cos_scenes'};
ROIs = cell(1, length(ROI_types)*2);
for r=1:length(ROI_types)
    ROIs{r*2 - 1} = strcat('lh_', ROI_types{r});
    ROIs{r*2} = strcat('rh_', ROI_types{r});
end


%% Track conditions for each subject
% ===== read master parfile ===== %
formats = {'dynacat', 'staticat', 'across'};
short_category_labels = {'words', 'dogs', 'bodies', 'hands', 'faces', 'cars', 'balls', 'scenes'};
category_labels = {'words-d', 'dogs-d', 'bodies-d', 'hands-d', 'faces-d', 'cars-d', 'balls-d', 'scenes-d'...
    'words-s', 'dogs-s', 'bodies-s', 'hands-s', 'faces-s', 'cars-s', 'balls-s', 'scenes-s'};
all_stimuli_cond_parfile_path = fullfile(ExpDir, 'code', 'experiment_code', 'DynaCat', 'dynacat_stimuli', 'phase_scrambled_set1', 'all_stimuli_conditions.par');
all_conds = {}; % record category for each condition number

fid = fopen(all_stimuli_cond_parfile_path, 'r');
data = textscan(fid, '%f %s %s %s %f %f %f', 'Delimiter', '\t');
fclose(fid);
for i=1:length(data{1})
    experiment = char(data{3}(i));
    category = char(data{4}(i));
    all_conds(end+1) = {sprintf('%s_%s',experiment, category)}; % e.g. dynacat_faces
end
for i = 49:72
    fprintf('%s\n', cell2str(data{2}(i)));
end
for i = 97:120
    fprintf('%s\n', cell2str(data{2}(i)));
end

% ===== read subject parfiles ===== %
cond_dict = containers.Map;
session_paths = {};
for s=1:length(subjects)
    subject = subjects{s};
    [session_path, list_runs, datatype, scan] = dynacat_staticat_sessions(subject);
    session_paths(end+1) = {session_path};
    cd(session_path);

    % ===== add individual stim parfiles ===== %
    parfiles = dir(fullfile(session_path, 'parfiles', 'individual_stim_parfiles_set1'));
    newscans_parfiles = {parfiles(3:end).folder}';
    newscans_parfiles = fullfile(newscans_parfiles, {parfiles(3:end).name}');

    % ===== find condition info from parfiles ===== %
    cond_nums = [];
    for i=1:length(newscans_parfiles)
        file_path = newscans_parfiles{i};
        fid = fopen(file_path, 'r');
        data = textscan(fid, '%s %f %s %f %f %f', 'Delimiter', '\t');
        fclose(fid);
        for ii=1:length(data{4})
            cond_num = data{2}(ii);
            if cond_num ~= 0 && cond_num ~= 999
                cond_nums = union(cond_nums, cond_num);
            end
        end
    end

    % convert condition numbers to words
    conditions = {};
    for i=1:length(cond_nums)
        conditions(end+1) = all_conds(cond_nums(i));
    end
    cond_dict(subject) = conditions;
end


%% Distinctiveness
distinctiveness_w = containers.Map; % within-format
distinctiveness_a = containers.Map; % across-format
distinctiveness_ds = containers.Map; % dynamic/static
distinctiveness_ai = containers.Map; % animate/inanimate
%load(fullfile(ResultsDir, 'scores.mat'));

for r=1:length(ROIs)
    ROI = ROIs{r};
    scores_all_subs_w = []; % grids of conditions across all subjects
    scores_all_subs_a = [];
    scores_all_subs_ds = [];
    scores_all_subs_ai = [];

    for s=1:length(subjects)
        subject = subjects{s};
        session_path = session_paths{s};
        if subjHasROI(session_path, ROI)
            conditions = cond_dict(subject);

            % within
            scores_w = calcDistinct_w(session_path, conditions, ROI);
            scores_all_subs_w = [scores_all_subs_w; scores_w];

            % across
            scores_a = calcDistinct_a(session_path, conditions, ROI);
            scores_all_subs_a = [scores_all_subs_a; scores_a];

            % dynamic/static
            scores_ds = calcDistinct_ds(session_path, conditions, ROI);
            scores_all_subs_ds = [scores_all_subs_ds; scores_ds];

            % animate/inanimate
            scores_ai = calcDistinct_ai(session_path, conditions, ROI);
            scores_all_subs_ai = [scores_all_subs_ai; scores_ai];

        end
    end
    distinctiveness_w(ROI) = scores_all_subs_w; 
    distinctiveness_a(ROI) = scores_all_subs_a;
    distinctiveness_ds(ROI) = scores_all_subs_ds;
    distinctiveness_ai(ROI) = scores_all_subs_ai;
end

fprintf('Done calculating distinctiveness\n')
save(fullfile(ResultsDir, 'scores.mat'), 'distinctiveness_w', 'distinctiveness_a', '-append');


%% Plotting
%load(fullfile(ResultsDir, 'scores.mat'));
mean_distinct = []; % save mean distinctiveness for each ROI in each format
mean_distinct_errors = [];

% individual ROI plots
for r=1:length(ROIs)
    ROI = ROIs{r};
    scores_all_subs_w = distinctiveness_w(ROI);
    scores_all_subs_a = distinctiveness_a(ROI);
    n_subs = size(scores_all_subs_w, 1);

    dynacat_means = mean(scores_all_subs_w(:, 1:8));
    staticat_means = mean(scores_all_subs_w(:,9:16));
    across_means = mean(scores_all_subs_a);
    data = [dynacat_means; staticat_means; across_means];
    mean_distinct = [mean_distinct mean(data, 2)];
    mean_distinct_errors = [mean_distinct_errors [std(data(1,:)); ...
        std(data(2,:)); std(data(3,:))]];

    hBar = bar(1:8, data);
    ylim([-0.05 0.4]);

    std_d = std(scores_all_subs_w(:, 1:8));
    std_s = std(scores_all_subs_w(:, 9:16));
    std_a = std(scores_all_subs_a);
    errors = [std_d; std_s; std_a];

    hold on
    nbars = size(data, 1);
    x = [];
    for n=1:nbars
        x = [x; hBar(n).XEndPoints];
    end
    errorbar(x, data, errors, 'k', 'linestyle', 'none');

    title_str = sprintf('%s (n = %d)', strrep(ROI, '_', ' '), n_subs);
    title(title_str);
    set(gca, 'xticklabel', short_category_labels);
    legend(formats, 'Location', 'Northwest');
    hold off

    filename = [ROI '_distinctiveness'];
    savefig(fullfile(FigDir, filename));
    saveas(gcf, fullfile(ImageDir, [filename '.png']));
    close

end

% graph all ROIs in each format
for i=1:length(formats)
    data = mean_distinct(i,:);
    data = [data(1:2:end); data(2:2:end)];

    errors = mean_distinct_errors(i,:);
    errors = [errors(1:2:end); errors(2:2:end)];

    ROI_types_labels = cell(1,length(ROI_types));
    for r=1:length(ROI_types)
        ROI_types_labels{r} = strrep(ROI_types{r}, '_', ' ');
    end

    hBar = bar(1:length(ROI_types), data);
    ylim([0 0.17]);

    hold on
    nbars = size(data, 1);
    x = [];
    for n=1:nbars
        x = [x; hBar(n).XEndPoints];
    end
    errorbar(x, data, errors, 'k', 'linestyle', 'none');

    title_str = sprintf('Distinctiveness across ROIs (%s)', formats{i});
    title(title_str);
    ax = hBar.Parent;
    set(ax, 'xtick', 1:length(ROI_types));
    set(gca, 'xticklabel', ROI_types_labels);
    set(gcf, 'position', [100, 100, 1700, 1000]);
    legend({'lh', 'rh'}, 'Location', 'Northwest');
    hold off

    filename = ['all_ROIs_' formats{i}];
    savefig(fullfile(FigDir, filename));
    saveas(gcf, fullfile(ImageDir, [filename '.png']));
    close
end


%% Plotting type comparisons
% Dynamic vs. Static
plotTypeComparison(ResultsDir, FigDir, ImageDir, ROI_types, 'format');

% Animate vs. Inanimate
plotTypeComparison(ResultsDir, FigDir, ImageDir, ROI_types, 'animacy');


%% Plotting everything
data = [];
errors = [];
for r=1:length(ROIs)
    ROI = ROIs{r};
    scores_a = mean(distinctiveness_a(ROI), 2);
    scores_ai = mean(distinctiveness_ai(ROI), 2);
    scores_ds = mean(distinctiveness_ds(ROI), 2);
    data = [data mean(scores_a) mean(scores_ai) mean(scores_ds)];
    errors = [errors std(scores_a) std(scores_ai) std(scores_ds)];
end

data = reshape(data, 6, [])';
data = data(:,[1 4 2 5 3 6]);
errors = reshape(errors, 6, [])';
errors = errors(:,[1 4 2 5 3 6]);

ROI_types_labels = cell(1,length(ROI_types));
for r=1:length(ROI_types)
    ROI_types_labels{r} = strrep(ROI_types{r}, '_', ' ');
end

hBar = bar(1:length(data), data);
ylim([-0.01 0.14]);
hBar(2).FaceColor = [0.5 0.75 0.93];
hBar(3).FaceColor = [0.85 0.35 0.1];
hBar(4).FaceColor = [1 0.625 0.5];
hBar(5).FaceColor = [0.5 0.18 0.6];
hBar(6).FaceColor = [0.75 0.6 0.8];

hold on
nbars = size(data, 2);
x = [];
for n=1:6
    x = [x; hBar(n).XEndPoints];
end
errorbar(x', data, errors, 'k', 'linestyle', 'none');

title('Distinctiveness across ROIs, all comparisons');
ax = hBar.Parent;
ylabel('Distinctiveness Score');
set(ax, 'xtick', 1:length(data));
set(gca, 'xticklabel', ROI_types_labels);
set(gcf, 'position', [100, 100, 2500, 1700]);
legend({'category', '', 'animate/inanimate', '', 'dynamic/static'}, 'Location', 'Northwest');
fontsize(28, 'points');
hold off

filename = 'the_big_plot';
savefig(fullfile(FigDir, filename));
saveas(gcf, fullfile(ImageDir, [filename '.png']));
close


%% Functions
function [has] = subjHasROI(session_path, ROI)
rsm_path = fullfile(session_path, 'results', 'func_rois_rsms', ...
    sprintf('%s_individual_stim_mean_rsm.mat', ROI));
has = isfile(rsm_path);
end

% returns whether current ptr is the end of a condition
function [at_end] = endOfSegment(ptr, rsm, conditions)
if ptr == length(rsm)
    at_end = 1;
else
    at_end = ~strcmp(conditions{ptr}, conditions{ptr+1});
end
end

% returns elements above the diagonal of a square matrix as an array
function [flattened] = aboveDiagonalFlat(M)
flattened = [];
for row=1:length(M)
    flattened = [flattened M(row, row+1:end)];
end
end

% find start of second format within subject rsm (staticat)
function [edge] = findFormatEdge(conditions)
edge = -1;
for c=1:length(conditions)
    if contains(conditions(c), 'staticat')
        edge = c;
        break;
    end
end
end

% distinctiveness scores for a subject within-format
function [scores] = calcDistinct_w(session_path, conditions, ROI)
scores = zeros(1, 16);
rsm_path = fullfile(session_path, 'results', 'func_rois_rsms', ...
    sprintf('%s_individual_stim_mean_rsm.mat', ROI));
rsm = load(rsm_path).rsm;

start = 1;
condition_idx = 1;
edge = findFormatEdge(conditions);
for ptr=1:length(rsm)
    if endOfSegment(ptr, rsm, conditions)
        % without cluster
        if start < edge % first/dynacat format
            without_cluster = rsm(start:ptr, [1:start-1, ptr+1:edge-1]);
        else % second/staticat format
            without_cluster = rsm(start:ptr, [edge:start-1, ptr+1:end]);
        end
        a = mean(mean(without_cluster));

        % within cluster
        within_cluster = aboveDiagonalFlat(rsm(start:ptr, start:ptr));
        b = mean(within_cluster);

        scores(condition_idx) = b - a; % 2 -> very distinct, -2 -> undistinct
        start = ptr + 1; % reset start idx
        condition_idx = condition_idx + 1; % next condition
    end
end
end

% distinctiveness scores for a subject across-format
function [scores] = calcDistinct_a(session_path, conditions, ROI)
    scores = zeros(1, 8);
    rsm_path = fullfile(session_path, 'results', 'func_rois_rsms', ...
        sprintf('%s_individual_stim_mean_rsm.mat', ROI));
    rsm = load(rsm_path).rsm;

    edge = findFormatEdge(conditions);
    start1 = 1; % first format
    start2 = edge; % second format
    condition_idx = 1;
    for ptr1=1:edge-1
        if endOfSegment(ptr1, rsm, conditions) % end of condition in first format
            % find end of condition in second format
            ptr2 = start2;
            while ptr2 < length(rsm) && strcmp(conditions(ptr2), conditions(ptr2+1))
                ptr2 = ptr2 + 1;
            end

            % without cluster
            without_cluster = rsm([start1:ptr1, start2:ptr2], ...
                [1:start1-1, ptr1+1:start2-1, ptr2+1:end]);
            a = mean(mean(without_cluster));

            % within cluster(s)
            within_cluster11 = aboveDiagonalFlat(rsm(start1:ptr1, start1:ptr1));
            within_cluster22 = aboveDiagonalFlat(rsm(start2:ptr2, start2:ptr2));
            within_cluster12 = rsm(start1:ptr1, start2:ptr2);
            b_sum = sum(within_cluster11) + sum(within_cluster22) + ...
                sum(sum(within_cluster12));
            b = b_sum / (length(within_cluster11) + ...
                length(within_cluster22) + numel(within_cluster12));

            scores(condition_idx) = b - a; % 2 -> very distinct, -2 -> undistinct
            start1 = ptr1 + 1; % reset starts
            start2 = ptr2 + 2;
            condition_idx = condition_idx + 1; % next condition
        end
    end
end

% distinctiveness scores between dynamic/static
function [scores] = calcDistinct_ds(session_path, conditions, ROI)
    scores = zeros(1, 2);
    rsm_path = fullfile(session_path, 'results', 'func_rois_rsms', ...
        sprintf('%s_individual_stim_mean_rsm.mat', ROI));
    rsm = load(rsm_path).rsm;

    edge = findFormatEdge(conditions);

    % without cluster (same for both formats)
    without_cluster = rsm(1:edge-1, edge:end);
    a = mean(mean(without_cluster));

    % within cluster for first format
    within_cluster1 = aboveDiagonalFlat(rsm(1:edge-1, 1:edge-1));
    b1 = mean(within_cluster1);
    scores(1) = b1 - a;

    % within cluster for second format
    within_cluster2 = aboveDiagonalFlat(rsm(edge:end, edge:end));
    b2 = mean(within_cluster2);
    scores(2) = b2 - a;
end

% distinctiveness scores between animacy/inanimacy
function [scores] = calcDistinct_ai(session_path, conditions, ROI)
    scores = zeros(1, 2);
    rsm_path = fullfile(session_path, 'results', 'func_rois_rsms', ...
        sprintf('%s_individual_stim_mean_rsm.mat', ROI));
    rsm = load(rsm_path).rsm;

    edge = findFormatEdge(conditions);
    start_dogs_d = -1; % start idx of animate categories in first format
    start_dogs_s = -1; % in second format
    start_cars_d = -1; % start idx of inanimate categories in first format
    start_cars_s = -1; % in second format
    idx = 1;
    while idx <= length(conditions) && start_cars_s == -1
        if start_dogs_d == -1 && strcmp(conditions(idx), 'dynacat_dogs')
            start_dogs_d = idx;
        end
        if start_dogs_s == -1 && strcmp(conditions(idx), 'staticat_dogs')
            start_dogs_s = idx;
        end
        if start_cars_d == -1 && strcmp(conditions(idx), 'dynacat_cars')
            start_cars_d = idx;
        end
        if start_cars_s == -1 && strcmp(conditions(idx), 'staticat_cars')
            start_cars_s = idx;
        end
        idx = idx + 1;
    end

    % without cluster (same for both formats)
    without_cluster = rsm([start_dogs_d:start_cars_d-1, ...
        start_dogs_s:start_cars_s-1], [1:start_dogs_d-1, ...
        start_cars_d:start_dogs_s-1, start_cars_s:end]);
    a = mean(mean(without_cluster));

    % within animate cluster
    dyn_to_dyn = aboveDiagonalFlat(rsm(start_dogs_d:start_cars_d-1, ...
        start_dogs_d:start_cars_d-1));
    sta_to_sta = aboveDiagonalFlat(rsm(start_dogs_s:start_cars_s-1, ...
        start_dogs_s:start_cars_s-1));
    dyn_to_sta = rsm(start_dogs_d:start_cars_d-1, start_dogs_s:start_cars_s-1);
    b1_sum = sum(dyn_to_dyn) + sum(sta_to_sta) + sum(sum(dyn_to_sta));
    b1 = b1_sum / (length(dyn_to_dyn) + length(dyn_to_sta) + numel(dyn_to_sta));
    scores(1) = b1 - a;

    % within inanimate cluster
    dyn_to_dyn = aboveDiagonalFlat(rsm(1:start_dogs_d-1, 1:start_dogs_d-1));
    dyn_to_dyn = [dyn_to_dyn aboveDiagonalFlat(rsm(start_cars_d:edge-1, ...
        start_cars_d:edge-1))];
    sta_to_sta = aboveDiagonalFlat(rsm(edge:start_dogs_s-1, edge:start_dogs_s-1));
    sta_to_sta = [sta_to_sta aboveDiagonalFlat(rsm(start_cars_s:end, ...
        start_cars_s:end))];
    dyn_to_sta_rest = rsm(start_cars_d:edge-1, start_cars_s:end);
    dyn_words_to_sta = rsm(1:start_dogs_d-1, start_cars_s:end);
    sta_words_to_dyn = rsm(edge:start_dogs_s-1, start_cars_d:edge-1);
    dyn_to_sta = [dyn_to_sta_rest(:)' dyn_words_to_sta(:)' sta_words_to_dyn(:)'];
    b2_sum = sum(dyn_to_dyn) + sum(sta_to_sta) + sum(dyn_to_sta);
    b2 = b2_sum / (length(dyn_to_dyn) + length(dyn_to_sta) + numel(dyn_to_sta));
    scores(2) = b2 - a;
end

% type – 'animacy','format'
function [] = plotTypeComparison(ResultsDir, FigDir, ImageDir, ROI_types, type)
load(fullfile(ResultsDir, 'scores.mat'));
if strcmp(type, 'animacy')
    vals = values(distinctiveness_ai);
    title_str = sprintf('Animacy distinctiveness across ROIs');
    legend_labels = {'lh (animacy)', 'rh (animacy)', 'lh (inanimacy)', 'rh (inanimacy)'};
elseif strcmp(type, 'format')
    vals = values(distinctiveness_ds);
    title_str = sprintf('Format distinctiveness across ROIs');
    legend_labels = {'lh (dynamic)', 'rh (dynamic)', 'lh (static)', 'rh (static)'};
else
    error('bad input');
end

data = [];
errors = [];
for i=1:length(vals)
    data = [data mean(vals{i})];
    errors = [errors std(vals{i})];
end
data = [data(1:4:end); data(2:4:end); data(3:4:end); ...
    data(4:4:end)];
errors = [errors(1:4:end); errors(2:4:end); errors(3:4:end); ...
    errors(4:4:end)];

ROI_types_labels = cell(1,length(ROI_types));
for r=1:length(ROI_types)
    ROI_types_labels{r} = strrep(ROI_types{r}, '_', ' ');
end

hBar = bar(1:length(data), data);
ylim([-0.01 0.12]);
hBar(2).FaceColor = [0.5 0.75 0.93];
hBar(3).FaceColor = [0.85 0.35 0.1];
hBar(4).FaceColor = [1 0.625 0.5];

hold on
nbars = size(data, 1);
x = [];
for n=1:nbars
    x = [x; hBar(n).XEndPoints];
end
errorbar(x, data, errors, 'k', 'linestyle', 'none');

title(title_str);
ax = hBar.Parent;
set(ax, 'xtick', 1:length(data));
set(gca, 'xticklabel', ROI_types_labels);
set(gcf, 'position', [100, 100, 1700, 1000]);
legend(legend_labels, 'Location', 'Northwest');
hold off

filename = [type '_distinctiveness'];
savefig(fullfile(FigDir, filename));
saveas(gcf, fullfile(ImageDir, [filename '.png']));
close
end