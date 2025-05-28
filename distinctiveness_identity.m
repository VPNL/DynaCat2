%% DynaCat + StatiCat: Identity
% Generate RSMs comparing the MVPs across identities of stimuli
%
% JC, February 2025


%% Setup
ExpDir = '/share/kalanit/biac2/kgs/projects/DynaCat/';
ImageDir = fullfile(ExpDir, 'results', 'anat_rois', 'scoring');
DataDir = fullfile(ExpDir, 'code', 'justin', 'rsms');

addpath(genpath(fullfile(ExpDir, 'code', 'scripts', 'rsms')));
subjects = {'AX', 'BR', 'CT', 'DO', 'HC', 'IK', 'KP', 'RU', 'RY', 'VL'};

ROI_types = {'IPS' , 'STS', 'LOTC', 'VTC', 'IF_DO'};
ROIs = cell(1, length(ROI_types)*2);
for r=1:length(ROI_types)
    ROIs{r*2 - 1} = strcat('lh_', ROI_types{r});
    ROIs{r*2} = strcat('rh_', ROI_types{r});
end

%% Parameters
nIdentities = 15;
identities = {'AX', 'BP', 'BR', 'BS', 'DO', 'EC', 'HO', 'IK', 'JC', 'JO', ...
    'JY', 'KP', 'KW', 'SE', 'SS'};

% store condition numbers where each identity is presented
condition_map = {[49 73 74 97 98] [50 51 75 76 99] [52 53 77 100 101] ...
    [54 78 102] [55 79 80 103 104] [56 57 81 105 106] [58 59 82 83 107 108] ...
    [60 61 84 109] [62 85 110] [63 64 86 87 88 111 112] [65 66 89 90 113 114] ...
    [67 68 91 92 93 115 116] [69 94 117 118] [70 71 95 119] [72 96 120]};

%% Calculate per-subject identity RSMs
for s = 1:length(subjects)
    subject = subjects{s};
    for r = 1:length(ROIs)
        ROI = ROIs{r};
        filename = [ROI '_rsm_all_cond.mat'];
        load(fullfile(DataDir, subject, filename));

        identity_rsm = nan(nIdentities, nIdentities);
        for i1 = 1:nIdentities
            for i2 = i1:nIdentities
                comparison = compScore(all_cond_rsm, condition_map{i1}, ...
                    condition_map{i2});
                identity_rsm(i1, i2) = comparison;
                identity_rsm(i2, i1) = comparison;
            end
        end

        % save out rsm
        rsm_filename = fullfile(DataDir, subject, [ROI '_identity_rsm.mat']);
        save(rsm_filename, 'identity_rsm');

        % plot as figure and save
        fig = plotRSM(identity_rsm, ROI, identities);
        rsm_png_filename = fullfile(DataDir, subject, [ROI '_identity_rsm.png']);
        saveas(fig, rsm_png_filename, 'png');
        close all;
    end
end

%% Generate across-subject identity RSM
for r = 1:length(ROIs)
    ROI = ROIs{r};
    big_rsm = nan(nIdentities, nIdentities, length(subjects));
    for s = 1:length(subjects)
        subject = subjects{s};
        rsm_filename = fullfile(DataDir, subject, [ROI '_identity_rsm.mat']);
        load(rsm_filename);
        big_rsm(:,:,s) = identity_rsm;
    end
    mean_rsm = nanmean(big_rsm, 3);

    % save out rsm
    rsm_filename = fullfile(DataDir, 'identity', [ROI '_identity_rsm.mat']);
    save(rsm_filename, 'mean_rsm');

    % plot as figure and save
    fig = plotRSM(identity_rsm, ROI, identities);
    rsm_png_filename = fullfile(DataDir, 'identity', [ROI '_identity_rsm.png']);
    saveas(fig, rsm_png_filename, 'png');
    close all;
end

%% Functions
% Compute the average correlation of identity MVPs given an RSM
function [score] = compScore(rsm, conditions1, conditions2)
sum_comps = 0;
total_comps = 0;
self_comparison = isequal(conditions1, conditions2);
if self_comparison
    i1_range = 1:length(conditions1) - 1;
else
    i1_range = 1:length(conditions1);
end

for i1 = i1_range
    % adjust range of inner loop
    if self_comparison
        i2_range = i1 + 1:length(conditions2);
    else
        i2_range = 1:length(conditions2);
    end

    for i2 = i2_range
        corr = rsm(conditions1(i1), conditions2(i2));
        if ~isnan(corr)
            sum_comps = sum_comps + corr;
            total_comps = total_comps + 1;
        end
    end
end

% account for edge case of no comparisons between identitites
if total_comps > 0
    score = sum_comps / total_comps;
else
    score = nan;
end

end

% Plot RSM
function [fig] = plotRSM(rsm, ROI, labels)
figure('color', [1 1 1], 'units', 'norm', 'position', [0.1 .1 .6 .6]);
imagesc(rsm, [-.7 .7]);
axis('image');
cmap = mrvColorMaps('coolhot');
colormap(cmap);
hcb = colorbar;
title(hcb, 'correlation');

plot_title = [ROI ': identity comparisons'];

title(plot_title, 'Interpreter','none', 'FontSize',24);
set(gca,'Xtick', 1:length(labels), 'XtickLabel', labels,'FontSize', 10);
set(gca,'Ytick', 1:length(labels), 'YtickLabel', labels, 'FontSize', 10);
set(gca, 'clim', [-0.35, 0.35]);
brighten(0.6)
fig = gcf;
end