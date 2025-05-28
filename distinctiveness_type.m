%% DynaCat + StatiCat: Distinctiveness Type
% Compute the distinctiveness of each type (without averaging)
%
% JC, April 2025

%% Setup
ExpDir = '/share/kalanit/biac2/kgs/projects/DynaCat/code/justin';
DataDir = fullfile(ExpDir, 'data');
FigDir = fullfile(ExpDir, 'figures');
RSMDir = fullfile(ExpDir, 'rsms');

addpath(genpath('/share/kalanit/biac2/kgs/projects/DynaCat/code/scripts/rsms'));
subjects = {'AX', 'BR', 'CT', 'DO', 'HC', 'IK', 'KP', 'RU', 'RY', 'VL'};

ROI_types = {'IPS', 'STS', 'LOTC', 'VTC', 'IF_DO'};
ROIs = cell(1, length(ROI_types)*2);
for r=1:length(ROI_types)
    ROIs{r*2 - 1} = strcat('lh_', ROI_types{r});
    ROIs{r*2} = strcat('rh_', ROI_types{r});
end

category_labels = ["words", "dogs", "bodies", "hands", "faces", "cars", ...
    "balls", "scenes"];
n_categories = length(category_labels);
n_conds = 384;
len_cond = 24;

hemis = {'lh', 'rh'};


%% Create Hypothesis Matrices
% Format Hypothesis
f_matrix = ones(n_conds, n_conds);
f_matrix(n_conds/2 + 1:end, 1:n_conds/2) = -1;
f_matrix(1:n_conds/2, n_conds/2 + 1:end) = -1;
for i = 1:n_conds
    f_matrix(i, i) = 0;
end

% Animacy Hypothesis
a_quadrant = ones(n_conds/2, n_conds/2);
a_quadrant(len_cond + 1:len_cond*5, [1:len_cond, len_cond*5 + 1:end]) = -1;
a_quadrant([1:len_cond, len_cond*5 + 1:end], len_cond + 1:len_cond*5) = -1;
a_matrix = repmat(a_quadrant, 2, 2);
for i = 1:n_conds
    a_matrix(i, i) = 0;
end

% Category Hypothesis
c_quadrant = -1*ones(n_conds/2, n_conds/2);
for i = 1:n_categories
    c_quadrant((i - 1)*len_cond + 1:i*len_cond, (i - 1)*len_cond + 1:i*len_cond) = 1;
end
c_matrix = repmat(c_quadrant, 2, 2);
for i = 1:n_conds
    c_matrix(i, i) = 0;
end


%% Compute distinctiveness
% store each distinctiveness score in a 3D matrix, indexed score type x ROI x
% subject
fprintf('Computing distinctiveness...\n');

n_types = 3;
scores_by_type = nan(n_types, length(ROIs), length(subjects));
for s = 1:length(subjects)
    subject = subjects{s};
    for r = 1:length(ROIs)
        ROI = ROIs{r};
        filename = [ROI '_rsm_all_cond.mat']; % load all_cond_rsm
        load(fullfile(RSMDir, subject, filename));

        % format, animacy, category
        scores_by_type(1, r, s) = distinctiveness(all_cond_rsm, f_matrix);
        scores_by_type(2, r, s) = distinctiveness(all_cond_rsm, a_matrix);
        scores_by_type(3, r, s) = distinctiveness(all_cond_rsm, c_matrix);
    end
end

cd(DataDir);
save('scores.mat', 'scores_by_type', '-append');
fprintf('Scores saved and calculated!\n\n');


%% Plotting
fprintf('Plotting distinctiveness per condition...\n');
cd(FigDir);
load(fullfile(DataDir, 'scores.mat'), 'scores_by_type');

for h = 1:length(hemis)
    hemi = hemis{h};
    hemi_scores = permute(scores_by_type(:,h:2:10,:), [2 1 3]);
    hemi_scores([1 4],:,:) = hemi_scores([4 1],:,:); % reorder ROIs for graph
    hemi_scores([2 3],:,:) = hemi_scores([3 2],:,:);

    means = mean(hemi_scores, 3);
    stds = std(hemi_scores, 0, 3);
    ROI_labels = {'VTC', 'LOTC', 'STS', 'IPS', 'VLPFC'};
    plotDistinctivenessByType(hemi, ROI_labels, means, stds);

    fprintf(sprintf('All done plotting %s!\n\n', hemi));
end


%% Statistics
fprintf('Running statistics on distinctiveness by condition...\n');
load(fullfile(DataDir, 'scores.mat'), 'scores_by_type');

for h = 1:length(hemis)
    hemi = hemis{h};

    % build the d-scores data structure
    n_scores = length(ROI_types)*3;
    d_scores = nan(length(subjects), n_scores); % subjects x d-scores
    for s = 1:length(subjects)
        for r = 1:length(ROI_types)
            d_scores(s, (r - 1)*3 + 1:r*3) = scores_by_type(:,r*2 - mod(h, 2),s)';
        end
    end

    col_names = strings([1, n_scores]);
    for i=1:n_scores
        col_names(i) = sprintf("d%d", i);
    end
    tbl = array2table(d_scores, VariableNames = col_names);

    % build the within-subject design structure
    condition = repmat(["format"; "animacy"; "category"], length(ROI_types), 1);
    ROI = repelem(string(ROI_types)', 3, 1);
    w2design = table(condition, ROI, VariableNames = ["condition", "ROI"]);

    % run the ANOVA
    rm = fitrm(tbl, "d1-d15~1", WithinDesign = w2design); % fit the RM model
    RT = ranova(rm, WithinModel = "condition * ROI");

    % clean data for figure
    pValue_clean = strings([length(RT.pValue) 1]);
    pValueGG_clean = strings([length(RT.pValueGG) 1]);
    F_clean = strings([length(RT.F) 1]);
    for i=1:length(RT.pValue)
        if mod(i, 2) == 0
            pValue_clean(i) = "";
            pValueGG_clean(i) = "";
            F_clean(i) = "";
        else
            pValue_clean(i) = num2str(RT.pValue(i, 1));
            pValueGG_clean(i) = num2str(RT.pValueGG(i, 1));
            F_clean(i) = num2str(RT.F(i, 1));
        end
    end
    RT = removevars(RT, ["SumSq", "MeanSq", "F", "pValue", "pValueGG", ...
        "pValueHF", "pValueLB"]);
    RT.F = F_clean;
    RT.pValue = pValue_clean;
    RT.pValueGG = pValueGG_clean;

    error_idx = contains(RT.Properties.RowNames, 'Error');
    RT(error_idx,:) = [];
    RT('(Intercept)',:) = [];
    row_names_clean = cell(size(RT, 1), 1);
    for i = 1:size(row_names_clean, 1)
        row_names_clean{i} = replace(RT.Properties.RowNames{i}, '(Intercept):', '');
    end
    RT.Properties.RowNames = row_names_clean;

    % plot
    filename = [upper(hemi) '_rm_anova'];
    title = [upper(hemi) ': Repeated Measures ANOVA'];
    fig = uifigure("Name", filename, "Position", [0 0 600 160]);
    uit = uitable(fig, "Data", RT(:,:), "Position", [20 20 550 110]);
    title_element = uilabel(fig, 'Text', title, 'Position', [225 130 350 20], ...
        'FontWeight', 'bold');

    % save as png
    SaveDir = fullfile(ExpDir, 'figures', 'anova');
    outfile = fullfile(SaveDir, [filename '.png']);
    exportapp(fig, outfile);
    close all
end

fprintf('All done!\n\n');

%% Functions
function vector = mat2vec(matrix)
vector = reshape(matrix, 1, []);
end

function res = distinctiveness(m, design)
total_b = 0;
total_a = 0;
n_elems_b = 0;
n_elems_a = 0;

for r = 1:size(m, 1)
    for c = 1:size(m, 2)
        if r ~= c && ~isnan(m(r, c))
            if design(r, c) == 1
                total_b = total_b + m(r, c);
                n_elems_b = n_elems_b + 1;
            else
                total_a = total_a + m(r, c);
                n_elems_a = n_elems_a + 1;
            end
        end
    end
end

b = total_b/n_elems_b;
a = total_a/n_elems_a;

res = b - a;
end

% dot product between two matrices + ignore NaNs + ignore diagonal +
% averaged over valid elements
function product = dot_product_avg(m1, m2)
total = 0;
n_elems = 0;
for r = 1:size(m1, 1)
    for c = 1:size(m1, 2)
        if r ~= c && ~isnan(m1(r, c)) && ~isnan(m2(r, c))
            total = total + (m1(r, c)*m2(r, c));
            n_elems = n_elems + 1;
        end
    end
end

product = total/n_elems;
end


% plots distinctivess by type (one hemisphere)
function [] = plotDistinctivenessByType(hemi, xlabs, data, errors)
hBar = bar(data, 'facecolor', 'flat');
hold on

ylim([0 0.1]);
yticks(0:0.02:0.1);

ngroups = size(data, 1);
nbars = size(data, 2);
groupwidth = min(0.8, nbars / (nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i - 1) * groupwidth / (2 * nbars);
    er = errorbar(x, data(:, i), errors(:, i), '.');
    er.Color = [0 0 0];
end

title(sprintf('%s: Distinctiveness by Condition\n', upper(hemi)), 'FontSize', 24);
ax = hBar.Parent;
set(ax, 'xtick', 1:length(data));
set(gca, 'xticklabel', xlabs, 'FontSize', 24);
set(gcf, 'position', [100, 100, 1700, 1000]);
xlabel(sprintf('\nROI'), 'FontWeight', 'bold', 'FontSize', 24);
ylabel(sprintf('Distinctiveness\n'), 'FontWeight', 'bold', 'FontSize', 24);
legend({'format', 'animacy', 'category'}, 'Location', 'northeast', 'FontSize', 24);

filename = [hemi '_distinctiveness_by_type.png'];
saveas(gcf, filename);
hold off
close
end
