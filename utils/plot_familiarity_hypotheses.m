%% DynaCat + Staticat: Plot Familiarity Hypotheses
%   
% BR, JC

%% Set parameters
FigDir = '/share/kalanit/biac2/kgs/projects/DynaCat/results/familiarity';
addpath(genpath('/share/kalanit/biac2/kgs/projects/DynaCat/code/scripts/rsms'));

n_conds = 12;
labels = {'bodies-d-unfamiliar', 'bodies-d-familiar', 'bodies-d-self', ...
    'bodies-s-unfamiliar', 'bodies-s-familiar', 'bodies-s-self', ...
    'faces-d-unfamiliar', 'faces-d-familiar', 'faces-d-self', ...
    'faces-s-unfamiliar', 'faces-s-familiar', 'faces-s-self'};

%% Plot hypothesess
hyp_titles = {'(a) Familiarity Hypothesis', '(b) Self Hypothesis', ...
    '(c) Category Hypothesis'};
heat = 1;
matrix = -heat * ones(n_conds, n_conds, 3);

% Familiarity Hypothesis
matrix([2 5 8 11], [2 5 8 11], 1) = heat;

% Self Hypothesis
matrix([3 6 9 12], [3 6 9 12], 2) = heat;

% Category Hypothesis
matrix(1:6, 1:6, 3) = heat;
matrix(7:12, 7:12, 3) = heat;

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
    set(gca, 'Xtick', 1:12, 'XtickLabel', labels, 'FontSize', 14);
    if i == 1 % y axis
        set(gca, 'Ytick', 1:12, 'YtickLabel', labels, 'FontSize', 14);
    else
        set(gca, 'Ytick', []);
    end
    if i == 3 % color bar
        hcb = colorbar; title(hcb, 'correlation');
        hcb.Layout.Tile = 'east';
        hcb.FontSize = 18;
    end
end

set(gcf, 'position', [100, 100, 1600, 600], 'color', 'white');

% save figure
filename = 'familiarity_hypotheses.png';
saveas(gcf, filename);
close all
