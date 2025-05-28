%% DynaCat + StatiCat: Distinctiveness (Repeated Measures ANOVA)
% Compute the distinctiveness of each category in
% each format, and then compile scores across subjects and run rANOVA.
%
% JC, April 2025

%% Setup
ExpDir = '/share/kalanit/biac2/kgs/projects/DynaCat/code/justin';
DataDir = fullfile(ExpDir, 'data');
SaveDir = fullfile(ExpDir, 'figures', 'anova');

addpath(genpath(fullfile(ExpDir, 'code', 'scripts', 'rsms')));
subjects = {'AX', 'BR', 'CT', 'DO', 'HC', 'IK', 'KP', 'RU', 'RY', 'VL'};

ROI_types = {'IPS', 'STS', 'LOTC', 'VTC', 'IFC'};
ROIs = cell(1, length(ROI_types)*2);
for r=1:length(ROI_types)
    ROIs{r*2 - 1} = strcat('lh_', ROI_types{r});
    ROIs{r*2} = strcat('rh_', ROI_types{r});
end

category_labels = ["words", "dogs", "bodies", "hands", "faces", "cars", ...
    "balls", "scenes"];
n_categories = length(category_labels);

%% Plotting distinctiveness per category
% load data
cd(DataDir);
load('scores.mat');

fprintf('Plotting distinctiveness per category...\n');
FigDir = fullfile(DataDir, 'figures', 'per_cat');
cd(FigDir);

for r = 1:length(ROI_types)
    comb_scores_d = zeros(n_categories, 2, length(subjects));
    comb_scores_s = zeros(n_categories, 2, length(subjects));

    lh_idx = 2*r - 1; % lh roi idx
    rh_idx = 2*r; % rh roi idx

    lh_scores_d = cat_scores_d(:,lh_idx,:);
    rh_scores_d = cat_scores_d(:,rh_idx,:);

    lh_scores_s = cat_scores_s(:,lh_idx,:);
    rh_scores_s = cat_scores_s(:,rh_idx,:);
    
    comb_scores_d(:,:,:) = [lh_scores_d rh_scores_d];
    comb_scores_s(:,:,:) = [lh_scores_s rh_scores_s];

    means = [mean(comb_scores_d, 3) mean(comb_scores_s, 3)];
    stds = [std(comb_scores_d, 0, 3) std(comb_scores_s, 0, 3)];
    plotDistinctivenessPerCat(ROI_types{r}, category_labels, means, stds);
end
fprintf('All done plotting!\n\n');


%% Three-way RM ANOVA
cd(DataDir);
load('scores.mat');

for r = 1:length(ROI_types)
    ROI = ROI_types{r}

    % build the d-scores data structure
    d_scores = nan(length(subjects), n_categories * 4); % subjects x d-scores
    d_scores_roi_d = cat_scores_d(:,2*r-1:2*r,:); % d-'code', 'justin', scores for roi
    d_scores_roi_s = cat_scores_s(:,2*r-1:2*r,:);
    for s=1:length(subjects)
        d_scores_roi_subj_d = d_scores_roi_d(:,:,s); % d-scores for roi x subj
        d_scores_roi_subj_s = d_scores_roi_s(:,:,s);
        d_scores(s,:) = [d_scores_roi_subj_d(:,1)' d_scores_roi_subj_d(:,2)' ...
            d_scores_roi_subj_s(:,1)' d_scores_roi_subj_s(:,2)'];
    end
    col_names = strings([1, 32]);
    for i=1:32
        col_names(i) = sprintf("d%d", i);
    end
    tbl = array2table(d_scores, VariableNames = col_names);

    % build the within-subject design structure
    format = repelem(["dynamic"; "static"], n_categories*2);
    hemi = repmat(repelem(["lh"; "rh"], n_categories), 2, 1);
    category = repmat(category_labels', 4, 1);
    w2design = table(format, hemi, category, VariableNames = ["format", ...
        "hemi", "category"]);

    % run the ANOVA
    rm = fitrm(tbl, "d1-d32~1", WithinDesign = w2design); % fit the RM model
    RT = ranova(rm, WithinModel = "format * hemi * category");

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
    row_names_clean = cell(7, 1);
    for i = 1:size(row_names_clean, 1)
        row_names_clean{i} = replace(RT.Properties.RowNames{i}, '(Intercept):', '');
    end
    RT.Properties.RowNames = row_names_clean;

    % plot
    if strcmp(ROI, 'IFC')
        ROI = 'VLPFC';
    end
    filename = [ROI '_rm_anova'];
    title = [ROI ': Repeated Measures ANOVA'];
    fig = uifigure("Name", filename, "Position", [0 0 600 250]);
    uit = uitable(fig, "Data", RT(:,:), "Position", [20 20 550 200]);
    title_element = uilabel(fig, 'Text', title, 'Position', [225 220 350 20], ...
        'FontWeight', 'bold');

    % save as png
    outfile = fullfile(SaveDir, [filename '.png']);
    exportapp(fig, outfile);
    close all
end


%% Functions
% plots distinctivess per category for a given ROI
function [] = plotDistinctivenessPerCat(ROI, xlabs, data, errors)
hBar = bar(data, 'facecolor', 'flat');
hold on

ylim([-0.03 0.24]);

colors = [0.2 0.647 0.9; 0 0.447 0.741; 1 0.425 0.198; 0.85 0.325 0.098];
for i = 1:8
    hBar(1).CData(i,:) = colors(1,:);
    hBar(2).CData(i,:) = colors(2,:);
    hBar(3).CData(i,:) = colors(3,:);
    hBar(4).CData(i,:) = colors(4,:);
end

ngroups = size(data, 1);
nbars = size(data, 2);
groupwidth = min(0.8, nbars / (nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i - 1) * groupwidth / (2 * nbars);
    er = errorbar(x, data(:, i), errors(:, i), '.');
    er.Color = [0 0 0];
end

title([replace(ROI, '_', '\_') ': Distinctiveness vs. Category'], 'FontSize', 16);
ax = hBar.Parent;
set(ax, 'xtick', 1:length(data));
set(gca, 'xticklabel', xlabs, 'FontSize', 16);
set(gcf, 'position', [100, 100, 1700, 1000]);
xlabel('\nCategory', 'FontWeight', 'bold', 'FontSize', 16);
ylabel('Distinctiveness\n', 'FontWeight', 'bold', 'FontSize', 16);
legend({'dynamic', '', 'static', ''}, 'Location', 'northwest');

filename = [ROI '_distinctiveness_per_cat_plot.png'];
saveas(gcf, filename);
hold off
close
end