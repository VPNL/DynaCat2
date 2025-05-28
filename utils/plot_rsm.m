%% DynaCat + Staticat: Plot RSMs
%   
% BR, JC

%% Set parameters
ExpDir = '/share/kalanit/biac2/kgs/projects/DynaCat/code/justin';
FigDir = fullfile(ExpDir, 'figures');

addpath(genpath('/share/kalanit/biac2/kgs/projects/DynaCat/code/scripts/rsms'));

ROI_types = {'IPS', 'STS', 'LOTC', 'VTC', 'IF_DO'};
ROIs = cell(1, length(ROI_types)*2);
for r = 1:length(ROI_types)
    ROIs{r} = strcat('lh_', ROI_types{r});
    ROIs{r + length(ROI_types)} = strcat('rh_', ROI_types{r});
end

n_conds = 384;
n_categories = 8;
len_cond = n_conds / (n_categories * 2);
category_labels = {'words-d', 'dogs-d', 'people-d', 'hands-d', 'faces-d', 'cars-d', 'balls-d', 'scenes-d'...
    'words-s', 'dogs-s', 'people-s', 'hands-s', 'faces-s', 'cars-s', 'balls-s', 'scenes-s'};


%% Plot all RSMs
cd(FigDir);

tiledlayout(2, 5, 'TileSpacing', 'tight');
tile_order = [4 3 2 1 5 9 8 7 6 10];
for i = 1:length(ROIs) % loop over ROIs
    r = tile_order(i);
    nexttile
    ROI = ROIs{r};
    ROIname = split(ROI, '_'); ROIname = ROIname(end);
    if strcmp(ROIname, 'DO')
        ROIname = 'VLPFC';
    end

    filename = [ROI '_rsm.mat'];
    load(fullfile(ExpDir, 'rsms', 'average', 'anatomical', filename)); % mean_rsm

    % plot rsm
    imagesc(mean_rsm, [-.7 .7]); axis('image');
    cmap = mrvColorMaps('coolhot'); colormap(cmap);
    set(gca, 'clim', [-1 1]);
    brighten(0.8)

    % figure labeling
    if i == 5 % color bar
        hcb = colorbar; title(hcb, 'correlation');
        hcb.Layout.Tile = 'east';
        hcb.FontSize = 18;
    end
    if i == 1 || i == 6 % y axis
        set(gca, 'Ytick', 1:24:max(n_conds), 'YtickLabel', category_labels, ...
            'FontSize', 14);
    else
        set(gca, 'Ytick', []);
    end
    if i <= 5 % ROI titles / x axis
        title(ROIname, 'FontSize', 30, 'FontWeight', 'normal');
        set(gca, 'Xtick', []);
    else
        set(gca,'Xtick', 1:24:max(n_conds), 'XtickLabel', category_labels, ...
            'FontSize', 14);
    end
    if i == 1
        ylabel('LH', 'FontSize', 30, 'FontWeight', 'bold');
    elseif i == 6
        ylabel('RH', 'FontSize', 30, 'FontWeight', 'bold');
    end
end

set(gcf, 'position', [100, 100, 1690, 700], 'color', 'white');

% save figure
filename = 'rsm_figure.png';
saveas(gcf, filename);
close all


%% Plot hypothesis
hyp_titles = {'(a) Format Hypothesis', '(b) Animacy Hypothesis', ...
    '(c) Category Hypothesis'};
heat = 1;
matrix = heat * ones(n_conds, n_conds, 3); % format, animacy, category

% Format Hypothesis
matrix(n_conds/2 + 1:end, 1:n_conds/2, 1) = -heat;
matrix(1:n_conds/2, n_conds/2 + 1:end, 1) = -heat;
for i = 1:n_conds
    matrix(i, i, 1) = 1;
end

% Animacy Hypothesis
a_quadrant = heat * ones(n_conds/2, n_conds/2);
a_quadrant(len_cond + 1:len_cond*5, [1:len_cond, len_cond*5 + 1:end]) = -heat;
a_quadrant([1:len_cond, len_cond*5 + 1:end], len_cond + 1:len_cond*5) = -heat;
matrix(:,:,2) = repmat(a_quadrant, 2, 2);
for i = 1:n_conds
    matrix(i, i, 2) = 1;
end

% Category Hypothesis
c_quadrant = -heat*ones(n_conds/2, n_conds/2);
for i = 1:n_categories
    c_quadrant((i - 1)*len_cond + 1:i*len_cond, (i - 1)*len_cond + 1:i*len_cond) = heat;
end
matrix(:,:,3) = repmat(c_quadrant, 2, 2);
for i = 1:n_conds
    matrix(i, i, 3) = 1;
end

% plot
cd(FigDir);
tiledlayout(1, 3, 'TileSpacing', 'tight');
for i = 1:size(matrix, 3) % loop over ROIs
    nexttile

    % plot rsm
    imagesc(matrix(:,:,i), [-.7 .7]); axis('image');
    cmap = mrvColorMaps('coolhot'); colormap(cmap);
    set(gca, 'clim', [-1 1]);
    brighten(0.6);

    % figure labeling
    ttl = title(sprintf('%s\n', hyp_titles{i}), 'FontSize', 36);
    set(gca, 'Xtick', []);
    set(gca, 'Xtick', 1:24:max(n_conds), 'XtickLabel', category_labels, ...
        'FontSize', 14);
    if i == 3 % color bar
        hcb = colorbar; title(hcb, 'correlation');
        hcb.Layout.Tile = 'east';
        hcb.FontSize = 18;
    end
    if i == 1 % y axis
        set(gca, 'Ytick', 1:24:max(n_conds), 'YtickLabel', category_labels, ...
            'FontSize', 14);
    else
        set(gca, 'Ytick', []);
    end
end

set(gcf, 'position', [100, 100, 1600, 600], 'color', 'white');

% save figure
filename = 'hypothesis_schematic.png';
saveas(gcf, filename);
close all
