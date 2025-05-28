%% DynaCat + StatiCat: Minigraphs
% Calculate + load distinctiveness of animacy + format across ROIs
%
% JC, April 2025

%% Setup
ExpDir = '/share/kalanit/biac2/kgs/projects/DynaCat/';
DataDir = fullfile(ExpDir, 'code', 'justin', 'rsms');
FigDir = fullfile(DataDir, 'figures', 'minigraphs');

addpath(genpath(fullfile(ExpDir, 'code', 'scripts', 'rsms')));
subjects = {'AX', 'BR', 'CT', 'DO', 'HC', 'IK', 'KP', 'RU', 'RY', 'VL'};

ROI_types = {'IPS', 'STS', 'LOTC', 'VTC', 'IF_DO'};
ROIs = cell(1, length(ROI_types)*2);
for r=1:length(ROI_types)
    ROIs{r*2 - 1} = strcat('lh_', ROI_types{r});
    ROIs{r*2} = strcat('rh_', ROI_types{r});
end

category_labels = {'words', 'dogs', 'bodies', 'hands', 'faces', 'cars', ...
    'balls', 'scenes'};
n_categories = length(category_labels);

%% Distinctiveness (animacy)
% store each distinctiveness in a 3D matrix, indexed animacy x ROI x
% subject
fprintf('Calculating distinctiveness (animacy)...\n');

animacy_scores_d = nan(2, length(ROIs), length(subjects));
animacy_scores_s = nan(2, length(ROIs), length(subjects));
animacy_scores_a = nan(2, length(ROIs), length(subjects));

for s = 1:length(subjects)
    subject = subjects{s};
    for r = 1:length(ROIs)
        ROI = ROIs{r};
        filename = [ROI '_rsm_all_cond.mat'];
        load(fullfile(DataDir, subject, filename));

        subj_scores = distinctiveness_animacy(all_cond_rsm);
        animacy_scores_d(:,r,s) = subj_scores(1,:)';
        animacy_scores_s(:,r,s) = subj_scores(2,:)';
        animacy_scores_a(:,r,s) = subj_scores(3,:)';
    end
end

cd(DataDir);
save('animacy_scores.mat', 'animacy_scores_d', 'animacy_scores_s', 'animacy_scores_a');
fprintf('Done calculating distinctiveness (animacy)!\n\n');


%% Distinctiveness (format)
% store each distinctiveness in a 3D matrix, indexed format x ROI x
% subject
fprintf('Calculating distinctiveness (format)...\n');

format_scores = nan(2, length(ROIs), length(subjects));

for s = 1:length(subjects)
    subject = subjects{s};
    for r = 1:length(ROIs)
        ROI = ROIs{r};
        filename = [ROI '_rsm_all_cond.mat'];
        load(fullfile(DataDir, subject, filename));

        subj_scores = distinctiveness_format(all_cond_rsm);
        format_scores(:,r,s) = subj_scores';
    end
end

cd(DataDir);
save('format_scores.mat', 'format_scores');
fprintf('Done calculating distinctiveness (format)!\n\n');


%% Plotting distinctiveness (animacy)
fprintf('Plotting distinctiveness (animacy)...\n');
cd(FigDir);

for r = 1:length(ROI_types)
    comb_scores_d = zeros(2, 2, length(subjects));
    comb_scores_s = zeros(2, 2, length(subjects));
    comb_scores_a = zeros(2, 2, length(subjects));

    lh_idx = 2*r - 1; % lh roi idx
    rh_idx = 2*r; % rh roi idx

    lh_scores_d = animacy_scores_d(:,lh_idx,:);
    rh_scores_d = animacy_scores_d(:,rh_idx,:);

    lh_scores_s = animacy_scores_s(:,lh_idx,:);
    rh_scores_s = animacy_scores_s(:,rh_idx,:);

    lh_scores_a = animacy_scores_a(:,lh_idx,:);
    rh_scores_a = animacy_scores_a(:,rh_idx,:);
    
    comb_scores_d(:,:,:) = [lh_scores_d rh_scores_d];
    comb_scores_s(:,:,:) = [lh_scores_s rh_scores_s];
    comb_scores_a(:,:,:) = [lh_scores_a rh_scores_a];

    means = [mean(comb_scores_d, 3) mean(comb_scores_s, 3) mean(comb_scores_a, 3)];
    stds = [std(comb_scores_d, 0, 3) std(comb_scores_s, 0, 3) std(comb_scores_a, 0, 3)];
    plotDistinctivenessAnimacy(ROI_types{r}, {'animate', 'inanimate'}, means, stds);
end
fprintf('All done plotting!\n\n');


%% Plotting distinctiveness (format)
fprintf('Plotting distinctiveness (format)...\n');
cd(FigDir);

for r = 1:length(ROI_types)
    roi_scores = permute(format_scores(:,[2*r-1 2*r],:), [2 1 3]);

    means = mean(roi_scores, 3);
    stds = std(roi_scores, 0, 3);
    plotDistinctivenessFormat(ROI_types{r}, {'dynamic', 'static'}, means, stds);
end
fprintf('All done plotting!\n\n');


%% Functions
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


% calculate distinctiveness of animacy + inanimacy
function [scores] = distinctiveness_animacy(rsm)
% across format
ani_d = aboveDiagonalFlat(rsm(25:120, 25:120));
ani_s = aboveDiagonalFlat(rsm(217:312, 217:312));
ani_a = flat(rsm(25:120, 217:312));

ina_words_d = aboveDiagonalFlat(rsm(1:24, 1:24));
ina_mixed_d = aboveDiagonalFlat(rsm(121:192, 121:192));
ina_words_to_mixed_d = flat(rsm(1:24, 121:192));
ina_words_s = aboveDiagonalFlat(rsm(193:216, 193:216));
ina_mixed_s = aboveDiagonalFlat(rsm(313:end, 313:end));
ina_words_to_mixed_s = flat(rsm(193:216, 313:end));
ina_a = flat(rsm([1:24, 121:192], [193:216, 313:end]));

without_d = flat(rsm(25:120, [1:24, 121:192]));
without_s = flat(rsm(217:312, [193:216, 313:end]));
without_a = [flat(rsm(25:120, [193:216, 313:end])) ...
    flat(rsm([1:24, 121:192], 193:216))];

scores_d = [(nanmean(ani_d) - nanmean(without_d)) ...
    (nanmean([ina_words_d ina_mixed_d ina_words_to_mixed_d]) - nanmean(without_d))];
scores_s = [(nanmean(ani_s) - nanmean(without_s)) ...
    (nanmean([ina_words_s ina_mixed_s ina_words_to_mixed_s]) - nanmean(without_s))];
scores_a = [(nanmean(ani_a) - nanmean(without_a)) (nanmean(ina_a) - nanmean(without_a))];

% animate | inanimate
scores = [scores_d; scores_s; scores_a];
end


% calculate distinctiveness of each format
function [scores] = distinctiveness_format(rsm)
dyn_to_dyn = nanmean(aboveDiagonalFlat(rsm(1:192, 1:192)));
sta_to_sta = nanmean(aboveDiagonalFlat(rsm(193:end, 193:end)));
dyn_to_sta = nanmean(flat(rsm(1:192, 193:end)));

distinctiveness_dyn = dyn_to_dyn - dyn_to_sta;
distinctiveness_sta = sta_to_sta - dyn_to_sta;

% dynamic | static
scores = [distinctiveness_dyn dcategory_labelsistinctiveness_sta];
end


% plot distinctivess (animacy) for a given ROI
function [] = plotDistinctivenessAnimacy(ROI, xlabs, data, errors)
hBar = bar(data, 'facecolor', 'flat');
hold on

ylim([-0.02 0.18]);
yticks([0:0.05:0.15]);

colors = [0.2 0.647 0.9; 0 0.447 0.741; 1 0.425 0.198; 0.85 0.325 0.098; ...
    1 0.85 0.5; 0.929 0.694 0.125];
for i = 1:2
    for c = 1:length(colors)
        hBar(c).CData(i,:) = colors(c,:);
    end
end

ngroups = size(data, 1);
nbars = size(data, 2);
groupwidth = min(0.8, nbars / (nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - granimacyoupwidth/2 + (2*i - 1) * groupwidth / (2 * nbars);
    er = errorbar(x, data(:, i), errors(:, i), '.');
    er.Color = [0 0 0];
end

if strcmp(ROI, 'IF_DO')
    ROI = 'IFC';
end
title(sprintf('%s: Distinctiveness vs. Animacy\n', ...
    replace(ROI, '_', '\_')), 'FontSize', 16);
ax = hBar.Parent;
set(ax, 'xtick', 1:length(data));
set(gca, 'xticklabel', xlabs, 'FontSize', 16);
set(gcf, 'position', [100, 100, 1700, 1000]);
xlabel(sprintf('\nAnimacy'), 'FontWeight', 'bold', 'FontSize', 16);
ylabel(sprintf('Distinctiveness\n'), 'FontWeight', 'bold', 'FontSize', 16);
legend({'dynamic', '', 'static', '', 'across', ''}, 'Location', 'northwest');

filename = [ROI '_distinctiveness_animacy_plot.png'];
saveas(gcf, filename);
hold off
close
end


% plot distinctivess (format18) for a given ROI
function [] = plotDistinctivenessFormat(ROI, xlabs, data, errors)
hBar = bar(data, 'facecolor', 'flat');
hold on

ylim([-0.005 0.05]);
yticks([0:0.01:0.04]);

colors = [0.2 0.647 0.9; 0 0.447 0.741; 1 0.425 0.198; 0.85 0.325 0.098];
hBar(1).CData(1,:) = colors(1,:);
hBar(2).CData(1,:) = colors(2,:);
hBar(1).CData(2,:) = colors(3,:);
hBar(2).CData(2,:) = colors(4,:);

ngroups = size(data, 1);
nbars = size(data, 2);
groupwidth = min(0.8, nbars / (nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i - 1) * groupwidth / (2 * nbars);
    er = errorbar(x, data(:, i), errors(:, i), '.');
    er.Color = [0 0 0];
end

if strcmp(ROI, 'IF_DO')
    ROI = 'IFC';
end
title(sprintf('%s: Distinctiveness vs. Format\n', ...
    replace(ROI, '_', '\_')), 'FontSize', 16);
ax = hBar.Parent;
set(ax, 'xtick', 1:length(data));
set(gca, 'xticklabel', xlabs, 'FontSize', 16);
set(gcf, 'position', [100, 100, 1700, 1000]);
xlabel(sprintf('\nFormat'), 'FontWeight', 'bold', 'FontSize', 16);
ylabel(sprintf('Distinctiveness\n'), 'FontWeight', 'bold', 'FontSize', 16);

filename = [ROI '_distinctiveness_format_plot.png'];
saveas(gcf, filename);
hold off
close
end
