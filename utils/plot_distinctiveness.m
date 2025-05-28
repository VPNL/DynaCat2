%% DynaCat + StatiCat: Plot Distinctiveness Graphs
% 1. Plot category distinctiveness across ROIs (grouped by format)
% 2. Plot condition distinctiveness per hemisphere
%
% JC, April 2025

%% Setup
ExpDir = '/share/kalanit/biac2/kgs/projects/DynaCat/code/justin';
DataDir = fullfile(ExpDir, 'data');
FigDir = fullfile(ExpDir, 'figures');

addpath(genpath('/share/kalanit/biac2/kgs/projects/DynaCat/code/scripts/rsms'));
subjects = {'AX' 'BR' 'CT' 'DO' 'HC' 'IK' 'KP' 'RU' 'RY' 'VL'};

ROI_types = {'IPS', 'STS', 'LOTC', 'VTC', 'IF_DO'};
ROIs = cell(1, length(ROI_types)*2);
for r=1:length(ROI_types)
    ROIs{r*2 - 1} = strcat('lh_', ROI_types{r});
    ROIs{r*2} = strcat('rh_', ROI_types{r});
end

category_labels = ["words", "dogs", "people", "hands", "faces", "cars", ...
    "balls", "scenes"];
n_categories = length(category_labels);
n_conds = 384;
len_cond = 24;

hemis = {'lh', 'rh'};


%% Plot category distinctiveness
load(fullfile(ExpDir, 'data', 'scores.mat'), 'cat_scores_d', 'cat_scores_s');
fprintf('Generating category distinctievness plot...\n');
SaveDir = fullfile(FigDir, 'category');
cd(SaveDir);

tiledlayout(1, 5, 'TileSpacing', 'tight');
tile_ordering = [4 3 2 1 5];
for i = 1:length(ROI_types)
    r = tile_ordering(i); % control order

    nexttile
    comb_scores_d = zeros(n_categories, 2, length(subjects));
    comb_scores_s = zeros(n_categories, 2, length(subjects));

    lh_idx = 2*r - 1; % lh roi idx
    rh_idx = 2*r; % rh roi idx
    
    comb_scores_d(:,:,:) = [cat_scores_d(:,lh_idx,:) cat_scores_d(:,rh_idx,:)];
    comb_scores_s(:,:,:) = [cat_scores_s(:,lh_idx,:) cat_scores_s(:,rh_idx,:)];

    means = [mean(comb_scores_d, 3) mean(comb_scores_s, 3)];
    stds = [std(comb_scores_d, 0, 3) std(comb_scores_s, 0, 3)];
    categoryDistinctivenessSubplot(ROI_types{r}, category_labels, means, stds);

    % y axis labels
    if i == 1
        ylabel(sprintf('Distinctiveness\n'), 'FontSize', 20, 'FontWeight', 'bold');
        set(gca, 'ytick', 0:0.1:0.3);
    else
        set(gca, 'ytick', []);
    end

    % x axis label
    if i == 3
        xlabel('Category', 'FontSize', 20, 'FontWeight', 'bold');
    end
    
end

% some figure elements
sgt = sgtitle(sprintf('Distinctiveness vs. Category\n'), 'FontWeight', 'bold');
sgt.FontSize = 24;
set(gcf, 'position', [100, 100, 2000, 400], 'color', 'white');
legend({'dynamic (lh)', 'dynamic (rh)', 'static (lh)', 'static (rh)'}, ...
    'FontSize', 12);

% save out
filename = 'category_distinctiveness.png';
saveas(gcf, filename);
hold off
close

fprintf('All done plotting!\n\n');


%% Plot condition distinctiveness
% load(fullfile(ExpDir, 'data', 'scores.mat'), 'scores_by_type');
fprintf('Generating condition distinctievness subplots...\n');
SaveDir = FigDir;
cd(SaveDir);

tiledlayout(2, 1, 'TileSpacing', 'tight');
for h = 1:length(hemis)
    nexttile
    set(gca, 'Visible', 'on');

    hemi = hemis{h};
    hemi_scores = permute(scores_by_type(:,h:2:10,:), [2 1 3]);
    hemi_scores([1 4],:,:) = hemi_scores([4 1],:,:); % reorder ROIs for graph
    hemi_scores([2 3],:,:) = hemi_scores([3 2],:,:);

    means = mean(hemi_scores, 3);
    stds = std(hemi_scores, 0, 3);
    if h == 1
        ROI_labels = {};
    else
        ROI_labels = {'VTC', 'LOTC', 'STS', 'IPS', 'VLPFC'};
    end
    plotDistinctivenessByType(ROI_labels, means, stds);

    if h == 1
        legend({'format', 'animacy', 'category'}, 'Location', 'northeast', 'FontSize', 16);
        ylabel(sprintf('LH'), 'FontWeight', 'bold', 'FontSize', 24);
    else
        xlabel(sprintf('\nROI'), 'FontWeight', 'bold', 'FontSize', 24);
        ylabel(sprintf('                                Distinctiveness\nRH'), ...
            'FontWeight', 'bold', 'FontSize', 24);
    end

end

sgtitle(sprintf('Distinctiveness by Condition\n'), 'FontWeight', 'bold', 'FontSize', 24);
set(gcf, 'position', [100, 100, 1000, 800], 'color', 'white');

% save out
filename = 'condition_distinctiveness.png';
saveas(gcf, filename);
hold off
close

fprintf('All done plotting!\n\n');


%% Functions
% create subplot of category distinctiveness
function [] = categoryDistinctivenessSubplot(ROI, xlabs, data, errors)
hBar = bar(data, 'facecolor', 'flat', 'EdgeColor', 'none');
hold on

ylim([-0.03 0.35]);

colors = [0.2 0.647 0.9; 0 0.447 0.741; 1 0.425 0.198; 0.85 0.325 0.098];
for i = 1:8
    for c = 1:length(colors)
        hBar(c).CData(i,:) = colors(c,:);
    end
end

ngroups = size(data, 1);
nbars = size(data, 2);
groupwidth = min(0.8, nbars / (nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i - 1) * groupwidth / (2 * nbars);
    errorbar(x, data(:, i), errors(:, i), '.', 'Color', 'black', 'LineWidth', 0.5);
end

box off
ax = hBar.Parent;
set(ax, 'xtick', 1:length(data));
set(gca, 'xticklabel', xlabs, 'FontSize', 16);
title(replace(ROI, 'IF_DO', 'VLPFC'), 'FontWeight', 'Normal', 'FontSize', 20);
end

% plots distinctivess by type (one hemisphere)
function [] = plotDistinctivenessByType(xlabs, data, errors)
hBar = bar(data, 'facecolor', 'flat');
hold on

ylim([0 0.14]);
yticks(0:0.04:0.12);

ngroups = size(data, 1);
nbars = size(data, 2);
groupwidth = min(0.8, nbars / (nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i - 1) * groupwidth / (2 * nbars);
    er = errorbar(x, data(:, i), errors(:, i), '.');
    er.Color = [0 0 0];
end

ax = hBar.Parent;
set(ax, 'xtick', 1:length(data));
set(gca, 'xticklabel', xlabs, 'FontSize', 24);
end
