%% DynaCat + StatiCat: Distinctiveness (ANOVA)
% Compute the distinctiveness of each category in
% each format, and then compile scores across subjects and run ANOVA.
%
% JC, March 2025

%% Setup
ExpDir = '/share/kalanit/biac2/kgs/projects/DynaCat/code/justin';
DataDir = fullfile(ExpDir, 'data');
FigDir = fullfile(ExpDir, 'figures');
ImageDir = fullfile(ExpDir, 'results', 'anat_rois', 'scoring');

addpath(genpath(fullfile(ExpDir, 'code', 'scripts', 'rsms')));
subjects = {'AX', 'BR', 'CT', 'DO', 'HC', 'IK', 'KP', 'RU', 'RY', 'VL'};

ROI_types = {'IPS', 'STS', 'LOTC', 'VTC', 'IF_DO'};
ROIs = cell(1, length(ROI_types)*2);
for r=1:length(ROI_types)
    ROIs{r*2 - 1} = strcat('lh_', ROI_types{r});
    ROIs{r*2} = strcat('rh_', ROI_types{r});
end

%category_labels_7 = {'words', 'dogs', 'hands', 'cars', 'balls', 'scenes', 'people'};
category_labels = {'words', 'dogs', 'bodies', 'hands', 'faces', 'cars', 'balls', ...
    'scenes'};
n_categories = 8;

%% Distinctiveness per category
% store each distinctiveness in a 3D matrix, indexed category x ROI x
% subject
% or in an array
fprintf('Calculating distinctiveness per category...\n');

%cat_scores_d = nan(n_categories, length(ROIs), length(subjects));
%cat_scores_s = nan(n_categories, length(ROIs), length(subjects));
cat_scores_a = nan(n_categories, length(ROIs), length(subjects));

for s = 1:length(subjects)
    subject = subjects{s};
    for r = 1:length(ROIs)
        ROI = ROIs{r};
        filename = [ROI '_rsm_all_cond.mat'];
        load(fullfile(DataDir, subject, filename));

        subj_scores = distinctiveness_per_category(all_cond_rsm);
        %scores_d = subj_scores(1,:);
        %scores_s = subj_scores(2,:);
        scores_a = subj_scores(3,:);

        cat_scores_a(:,r,s) = scores_a';
    end
end

cd(DataDir);
save('scores.mat', 'cat_scores_a', '-append');
fprintf('Done calculating distinctiveness per category!\n\n');


%% Four-way ANOVA
% hemi x ROI x category x format
cd(DataDir);
load('scores.mat');
nROIs = length(ROIs);
nSubjects = length(subjects);
g_category = repmat(category_labels, 1, nSubjects * nROIs * 2);
d_string = repmat({'d'}, 1, nCategories);
s_string = repmat({'s'}, 1, nCategories);
g_format = repmat([d_string, s_string], 1, nSubjects * nROIs);
lh_string = repmat({'lh'}, 1, nCategories * 2);
rh_string = repmat({'rh'}, 1, nCategories * 2);
g_hemi = repmat([lh_string, rh_string], 1, nSubjects * nROIs / 2);
ROIs_string = repelem(ROI_types, nCategories * 4);
g_ROI = repmat(ROIs_string, 1, nSubjects);

%p = anovan(cat_scores, {g_hemi, g_ROI, g_category, g_format})
p = anovan(cat_scores, {g_hemi', g_ROI', g_category', g_format'}, 'model', ...
    'interaction', 'varnames', {'hemi', 'ROI', 'category', 'format'})
%p = anovan(cat_scores, {g_hemi, g_ROI, g_category, g_format}, 3, 3, ...
%   {'hemi', 'ROI', 'category', 'format'});

%% Three-way ANOVAs
for r=1:length(ROI_types)
    ROI = ROI_types{r}

    % build the data structure
    cat_scores_r = [];
    g_hemi_r = [];
    g_category_r = [];
    g_format_r = [];
    for i=1:length(cat_scores)
        if strcmp(ROI, g_ROI{i})
            cat_scores_r = [cat_scores_r; cat_scores(i)];
            g_hemi_r = [g_hemi_r; g_hemi(i)];
            g_category_r = [g_category_r; g_category(i)];
            g_format_r = [g_format_r; g_format(i)];
        end
    end

    % run the ANOVA
    p = anovan(cat_scores_r, {g_hemi_r, g_category_r, g_format_r}, 'model', ...
        'interaction', 'varnames', {'hemi', 'category', 'format'});
    close
end

%% Functions
% returns elements above the diagonal of a square matrix as an array
function [flattened] = above_diagonal_flat(M)
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
cat_scores_a = zeros(1, 8);
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
        start_dyn + 324:192]);
    a_a = nanmean([flat(dyn_to_sta_other) flat(sta_to_dyn_other)]);
    a_d = nanmean(flat(dyn_to_other));
    a_s = nanmean(flat(sta_to_other));

    % within cluster
    dyn_cluster = above_diagonal_flat(rsm(start_dyn:start_dyn + 23, start_dyn:start_dyn + 23));
    sta_cluster = above_diagonal_flat(rsm(start_sta:start_sta + 23, start_sta:start_sta + 23));
    across_format_within_cat = rsm(start_dyn:start_dyn + 23, start_sta:start_sta + 23);
    b_a = nanmean(flat(across_format_within_cat));
    b_d = nanmean(dyn_cluster);
    b_s = nanmean(sta_cluster);

    cat_scores_a(condition) = (b_a - a_a);
    cat_scores_d(condition) = (b_d - a_d);
    cat_scores_s(condition) = (b_s - a_s);
end

scores = [cat_scores_d; cat_scores_s; cat_scores_a];
end


% calculates the distinctiveness scores of 7 categories in 2 formats, where
% the 7th category is person (bodies + faces)
%{
function [scores] = distinctiveness_per_category(rsm)
cat_scores_d = zeros(1, 7);
cat_scores_s = zeros(1, 7);

% for single, ungrouped categories record their indices in the original rsm
single_index_map = [1 2 4 6 7 8];
for cnd=1:6
    idx = single_index_map(cnd);
    start_d = ((idx - 1) * 24) + 1;
    start_s = start_d + 192;

    % without clustersize(colors,2)
    dyn_to_other = rsm(start_d:start_d + 23, [1:start_d - 1, ...
        start_d + 24:192]);
    sta_to_other = rsm(start_s:start_s + 23, [193:start_s - 1, ...
        start_s + 24:end]);
    a_d = nanmean(flat(dyn_to_other));
    a_s = nanmean(flat(sta_to_other));

    % within cluster
    dyn_cluster = above_diagonal_flat(rsm(start_d:start_d + 23, start_d:start_d + 23));
    sta_cluster = above_diagonal_flat(rsm(start_s:start_s + 23, start_s:start_s + 23));
    b_d = nanmean(dyn_cluster);
    b_s = nanmean(sta_cluster);

    cat_scores_d(cnd) = b_d - a_d;
    cat_scores_s(cnd) = b_s - a_s;
end

% calculate d-score for person cat
body_start_d = 49;
body_start_s = body_start_d + 192;
face_start_d = 97;
face_start_s = face_start_d + 192;

% without cluster
dyn_to_other = rsm([body_start_d:body_start_d + 23, face_start_d:face_start_d + 23], ...
    [1:body_start_d - 1, body_start_d + 24:face_start_d - 1, face_start_d + 24:192]);
sta_to_other = rsm([body_start_s:body_start_s + 23, face_start_s:face_start_s + 23], ...
    [1:body_start_s - 1, body_start_s + 24:face_start_s - 1, face_start_s + 24:end]);
a_d = nanmean(flat(dyn_to_other));
a_s = nanmean(flat(sta_to_other));

% within cluster
body_to_body_d = above_diagonal_flat(rsm(body_start_d:body_start_d + 23, ...
    body_start_d:body_3start_d + 23));
face_to_face_d = above_diagonal_flat(rsm(face_start_d:face_start_d + 23, ...
    face_start_d:face_start_dlength(colors) + 23));
body_to_face_d = flat(rsm(body_start_d:body_start_d + 23, ...
    face_start_d:face_start_d + 23));
dyn_cluster = [body_to_body_d face_to_face_d body_to_face_d];
b_d = nanmean(dyn_cluster);

body_to_body_s = above_diagonal_flat(rsm(body_start_s:body_start_s + 23, ...
    body_start_s:body_start_s + 23));
face_to_face_s = above_diagonal_flat(rsm(face_start_s:face_start_s + 23, ...
    face_start_s:face_start_s + 23));
body_to_face_s = flat(rsm(body_start_s:body_start_s + 23, ...
    face_start_s:face_start_s + 23));
sta_cluster = [body_to_body_s face_to_face_s body_to_face_s];
b_s = nanmean(sta_cluster);

cat_scores_d(7) = b_d - a_d;
cat_scores_s(7) = b_s - a_s;

scores = [cat_scores_d; cat_scores_s];
end
%}
