%% DynaCat + Staticat: Plot MDS
%   -Plot MDS for each ROI using pre-computed values
%
% BR, JC

%% Set parameters
DynaCat = '/share/kalanit/biac2/kgs/projects/DynaCat/';
DataDir = fullfile(DynaCat, 'results', 'individual_stimuli_analysis', 'across_subjects');
ExpDir = fullfile(DynaCat, 'code', 'justin');
FigDir = fullfile(ExpDir, 'figures');
addpath(genpath('/share/kalanit/biac2/kgs/projects/DynaCat/code/scripts/rsms/'))

ROI_types = {'IPS', 'STS', 'LOTC', 'VTC', 'IF_DO'};
ROIs = cell(1, length(ROI_types)*2);
for r=1:length(ROI_types)
    ROIs{r*2 - 1} = strcat('lh_', ROI_types{r});
    ROIs{r*2} = strcat('rh_', ROI_types{r});
end

stimuli_icon_folder = '/share/kalanit/biac2/kgs/projects/DynaCat/stimuli/stimuli_icons/stimuli_icons_lower_res/dynacat_staticat_colored_icons';
stimuli_icon_label = 'exp_colors';
stimuli_icon_size = 0.015;
categeory_color_point_size = 17;

num_all_conds = 384;
parfile_folders = 'individual_stim_parfiles_set1';

category_labels = {'words-d', 'dogs-d', 'bodies-d', 'hands-d', 'faces-d', 'cars-d', 'balls-d', 'scenes-d'...
    'words-s', 'dogs-s', 'bodies-s', 'hands-s', 'faces-s', 'cars-s', 'balls-s', 'scenes-s'};

% load stimuli condition info (for group mds plot labels)
all_stimuli_cond_parfile_path = fullfile(DynaCat, 'code', 'experiment_code', 'DynaCat', 'dynacat_stimuli', 'phase_scrambled_set1', 'all_stimuli_conditions.par');
all_cond_parfile = {};
data = readtable(all_stimuli_cond_parfile_path, 'FileType', 'delimitedtext');
for i = 1:height(data)
    cond_num = data{i, 1};
    cond = data{i, 2};
    experiment = data{i, 3};
    category = data{i, 4};
    color_1 = data{i, 5};
    color_2 = data{i, 6};
    color_3 = data{i, 7};
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

%% Plot stimulus icons
fprintf('Generating MDS plot (stimuli icons)...\n');
SaveDir = fullfile(FigDir, 'mds');
cd(SaveDir);

tiledlayout(5, 2);
tile_ordering = [7 8 5 6 3 4 1 2 9 10];
for idx = 1:length(ROIs)
    nexttile
    r = tile_ordering(idx);
    ROIname = ROIs{r};
    ROIlabel = split(ROIname, '_');
    ROIlabel = ROIlabel(end);
    if strcmp(ROIlabel, 'DO')
        ROIlabel = 'VLPFC';
    end
    fprintf('Plotting %s\n', ROIname);

    % get mean MDS across-subjects per ROI
    load(fullfile(DataDir, [ROIname '_mean_distance_matrix.mat'])); % mean_distance_matrix
    Y = cmdscale(mean_distance_matrix, 2);
    if idx == 1
        prosM = Y;
    else
        [d, Y, t] = procrustes(prosM, Y);
    end
    
    % plot each condition
    stimuli_icon_size_roi = 0.03;
    hold on 
    for i = 1:numel(all_cond_labels)
        if strcmp(ROIlabel, 'VTC') || strcmp(ROIname, 'rh_IPS')
            Y(i,1) = -Y(i,1);
        end
        scatter(Y(i,1), Y(i,2), 'Marker', 'none');

        [stim_icon, map, alphachannel] = imread(fullfile(stimuli_icon_folder, [all_cond_labels{i} '.png']));
        stim_icon = flipud(stim_icon); 
        alphachannel = flipud(alphachannel); 

        dx = stimuli_icon_size_roi; dy = stimuli_icon_size_roi;
        xmin = Y(i, 1) - dx; xmax = Y(i, 1) + dx;
        ymin = Y(i, 2) - dy; ymax = Y(i, 2) + dy;
        h = image([xmin xmax], [ymin ymax], stim_icon, 'AlphaData', alphachannel);
    end

    xlim([-0.6 0.6]); ylim([-0.6 0.6]);

    % figure labels
    if idx == 1
        title('LH', 'FontSize', 16);
    elseif idx == 2
        title('RH', 'FontSize', 16);
    end
    if mod(idx, 2) == 0
        set(gca, 'Ytick', []);
    else
        ylabel(ROIlabel, 'FontSize', 16);
    end
    if idx <= 8
        set(gca, 'Xtick', []);
    end
    hold off

    clear mean_distance_matrix
end

sgtitle(sprintf('MDS Plots (stimulus icons)\n'), 'FontWeight', 'bold', 'FontSize', 18);
set(gcf, 'position', [100, 100, 400, 800], 'color', 'white');

% save out
filename = 'mds_plot_stimulus_icons.png';
saveas(gcf, filename);
hold off
close

fprintf('All done!\n\n');


%% Plot category icons
fprintf('Generating MDS plot (category icons)...\n');
SaveDir = fullfile(FigDir, 'mds');
cd(SaveDir);

tiledlayout(5, 2);
tile_ordering = [7 8 5 6 3 4 1 2 9 10];
for idx = 1:length(ROIs)
    nexttile
    r = tile_ordering(idx);
    ROIname = ROIs{r};
    ROIlabel = split(ROIname, '_');
    ROIlabel = ROIlabel(end);
    if strcmp(ROIlabel, 'DO')
        ROIlabel = 'VLPFC';
    end
    fprintf('Plotting %s\n', ROIname);

    % get mean MDS across-subjects per ROI
    load(fullfile(DataDir, [ROIname '_mean_distance_matrix.mat'])); % mean_distance_matrix
    [Y, eigvals] = cmdscale(mean_distance_matrix);
    
    % make category colors mds
    categeory_color_point_size_roi = 10;
    hold on 
    for i = 1:numel(all_cond_labels)
        if strcmp(ROIlabel, 'VTC') || strcmp(ROIname, 'rh_IPS')
            Y(i,1) = -Y(i,1);
        end
        if i < 193
            scatter(Y(i,1), Y(i,2), categeory_color_point_size_roi, 'Marker', all_cond_point_shapes{i}, 'MarkerFaceColor', all_cond_label_colors{i}, 'MarkerEdgeColor', all_cond_label_colors{i});
        else
            scatter(Y(i,1), Y(i,2), categeory_color_point_size_roi, 'Marker', all_cond_point_shapes{i}, 'MarkerEdgeColor', all_cond_label_colors{i});
        end
    end
    xlim([-0.6 0.6]); ylim([-0.6 0.6]);

    xlim([-0.6 0.6]); ylim([-0.6 0.6]);

    % figure labels
    if idx == 1
        title('LH', 'FontSize', 16);
    elseif idx == 2
        title('RH', 'FontSize', 16);
    end
    if mod(idx, 2) == 0
        set(gca, 'Ytick', []);
    else
        ylabel(ROIlabel, 'FontSize', 16);
    end
    if idx <= 8
        set(gca, 'Xtick', []);
    end
    hold off

    clear mean_distance_matrix
end

sgtitle(sprintf('MDS Plots (category icons)\n'), 'FontWeight', 'bold', 'FontSize', 18);
set(gcf, 'position', [100, 100, 400, 800], 'color', 'white');

% save out
filename = 'mds_plot_category_icons.png';
saveas(gcf, filename);
hold off
close

fprintf('All done!\n\n');