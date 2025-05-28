%% DynaCat + StatiCat: Distinctiveness
% Compute within-format and across-format distinctiveness
%
% JC, November 2024


%% Setup
ExpDir = '/share/kalanit/biac2/kgs/projects/DynaCat/';
ImageDir = fullfile(ExpDir, 'results', 'anat_rois', 'scoring');
addpath(genpath(fullfile(ExpDir, 'code', 'scripts', 'rsms')));
subjects = {'AX', 'BR', 'CT', 'DO', 'HC', 'IK', 'KP', 'RU', 'RY', 'VL'};
formats = {'category', 'animacy', 'motion'}; % 'identity'

ROI_types = {'IPS' , 'STS', 'LOTC', 'VTC'};
ROIs = cell(1, length(ROI_types)*2);
for r=1:length(ROI_types)
    ROIs{r*2 - 1} = strcat('lh_', ROI_types{r});
    ROIs{r*2} = strcat('rh_', ROI_types{r});
end

prop_voxels = 0.1:0.1:1;
distinctiveness = containers.Map();
mean_ve = containers.Map();

% load data if calculated
load('data.mat');

%% GLM parameters
eventsPerBlock = 4; % length of trial is 4s and the TR is 1s
ampType = 'zscore';
plot_type = 'sym_means';

num_all_conds = 384;
parfile_folders = 'individual_stim_parfiles_set1';

category_labels = {'words-d', 'dogs-d', 'bodies-d', 'hands-d', 'faces-d', 'cars-d', 'balls-d', 'scenes-d'...
    'words-s', 'dogs-s', 'bodies-s', 'hands-s', 'faces-s', 'cars-s', 'balls-s', 'scenes-s'};

%% Compute subject-level MVPs + RSMs
for r = 1:length(ROIs)

    for v = 1:length(prop_voxels)
        ROIname = ROIs{r};
        prop = prop_voxels(v);
        ve_list = [];
        fprintf('%s at %.2f prop of voxels\n\n', ROIname, prop);

        % initialize across-subject rsm
        across_subject_rsm = nan(num_all_conds, num_all_conds, length(subjects));

        for s = 1:length(subjects)
            % ===== get subject session info =====
            subject = subjects{s};
            [session_path, list_runs, datatype, scan] = dynacat_staticat_sessions(subject);
            cd(session_path)

            % ===== add individual stim parfiles =====
            parfiles = dir(fullfile(session_path, 'parfiles', parfile_folders));
            newscans_parfiles = {parfiles(3:end).folder}';
            newscans_parfiles = fullfile(newscans_parfiles, {parfiles(3:end).name}');

            % ===== find condition info from parfiles =====
            parfile_data = {};
            for i=1:length(newscans_parfiles)
                file_path = newscans_parfiles{i};
                fid = fopen(file_path, 'r');
                data = textscan(fid, '%s %f %s %f %f %f', 'Delimiter', '\t');
                fclose(fid);
                for ii=1:length(data{4})
                    color_1 = data{4}(ii);
                    color_2 = data{5}(ii);
                    color_3 = data{6}(ii);
                    color_code = [color_1 color_2 color_3];
                    onset = str2double(data{1}(ii));
                    cond_num = data{2}(ii);
                    cond = data{3}(ii);
                    line_data = [i onset cond_num cond color_code];
                    parfile_data = [parfile_data; line_data];
                end
            end
            sorted_parfile_data = sortrows(parfile_data, 3); % sort parfile data based on condition number

            % remove baseline and task
            cond_nums = cell2mat(sorted_parfile_data(:, 3));
            rowsToKeep = (cond_nums ~= 0) & (cond_nums ~= 999);
            filtered_data = sorted_parfile_data(rowsToKeep, :);

            % remove some repeating condition labels that snuck their way in
            unique_items = {};
            rows_to_keep = true(size(filtered_data, 1), 1);
            for i = 1:size(filtered_data, 1)
                item = filtered_data{i, 4};
                if ismember(item, unique_items)
                    rows_to_keep(i) = false;
                else
                    unique_items = [unique_items; item];
                end
            end

            filtered_parfile_data = filtered_data(rows_to_keep, :);
            parfile_cond_nums = cell2mat(filtered_parfile_data(:, 3));
            numconds = length(unique(parfile_cond_nums));

            % ===== plot rsm =====
            % copy ROI over from mrVistaROI folder if nesecary
            ROI_filepath = fullfile(session_path, '3DAnatomy', 'ROIs', [ROIname '.mat']);
            ROI_filepath_in_mrVistaROI_dir = fullfile(session_path, '3DAnatomy', 'mrVistaROIs', [ROIname '.mat']);
            if ~exist(ROI_filepath, 'file')
                copyfile(ROI_filepath_in_mrVistaROI_dir, ROI_filepath);
            end

            % set vANATOMY
            hg = initHiddenGray(datatype, scan, ROIname);

            % add individual stim parfiles
            runs = list_runs;
            hg = er_assignParfilesToScans(hg, runs, newscans_parfiles);

            % initialize multivoxel pattern
            mv = mv_init(hg,ROIname,runs,datatype);
            mvs = {mv};

            % set parameters
            params = mvs{1}.params;
            params.eventAnalysis = 1;
            params.detrend = 1; % high-pass filter (remove the low frequency trend)
            params.detrendFrames = 20;
            params.inhomoCorrection = 1; % divide by mean (transforms from raw scanner units to % signal)
            params.temporalNormalization = 0; % no matching 1st temporal frame
            params.ampType = 'betas';
            params.glmHRF = 3; % SPM hrf
            params.eventsPerBlock = eventsPerBlock;
            params.lowPassFilter = 0; % no temporal low pass filter of the data

            mv.params = params;

            % apply GLM to multivoxel pattern
            mv = mv_applyGlm(mv);

            % get betas for each run and voxel
            amps = mv_amps(mv);

            % removes task condition from amps
            amps = amps(:, 1:numconds);

            % values to plot
            z_values = amps;

            if strcmp(ampType, 'zscore')
                % subtract mean beta from each voxel
                mean_amps=mean(amps,2); % calculate the mean for each voxel
                meanmat=mean_amps*ones(1,numconds);
                subbetas=amps-meanmat;
                % normalize by the residual variance of the GLM in each voxel
                % effectively after this step you have z-scoed your pattern
                residualV = mv.glm.residual; % make this into functions
                dof = mv.glm.dof;
                residualvar = sum(residualV.^2)/dof; % calculate mean sqaure error
                resd = sqrt(residualvar);
                resd = resd';
                resd_mat = repmat(resd, [1 numconds]); %turn into matrix
                z_values = subbetas./resd_mat; % divide each beta by standard deviation of residual error
            end

            % titrate
            nvoxels = round(prop * length(mv.glm.varianceExplained));
            [B, voxels_idx] = maxk(mv.glm.varianceExplained, nvoxels);
            z_values = z_values(voxels_idx',:);
            ve_list = [ve_list mean(B)];

            % compute crosscorrelation matrices
            rsm = zeros(numconds);
            all_cond_rsm = NaN(num_all_conds, num_all_conds);
            for i=1:numconds % calculate correlation between each condition
                cond_num_x = cell2mat(filtered_parfile_data(i, 3));
                for j=1:numconds
                    c = corrcoef(z_values(:,i),z_values(:,j));
                    rsm(i, j) = c(1, 2); % subject rsm

                    cond_num_y = cell2mat(filtered_parfile_data(j, 3));
                    all_cond_rsm(cond_num_x, cond_num_y) = rsm(i, j); % all cond subject rsm
                end
            end

            across_subject_rsm(:, :, s) = all_cond_rsm; % add this subject to across subject rsm

            % clear, close everything to save memory
            close all
            clear roi_rsm
            clear rsm
            clear all_cond_rsm
            clear mv
            clear z_values
            clear amps
            clear allvals
            clear cmap
            clear mean_amps
            clear meanmat
            clear mvs
            clear resd
            clear resd_mat
            clear residualV
            clear residualvar
            clear subbetas
        end

        % get mean across subjects rsms per roi
        mean_rsm = nanmean(across_subject_rsm, 3);

        % calculate distinctiveness
        d_score(1) = distinctiveness_category(mean_rsm);
        d_score(2) = distinctiveness_animacy(mean_rsm);
        d_score(3) = distinctiveness_motion(mean_rsm);

        % define maps
        for f = 1:length(formats) % distinctiveness across prop ve
            format = formats{f};
            if ~ismember(format, distinctiveness.keys)
                distinctiveness(format) = containers.Map;
            end

            format_map = distinctiveness(format);
            if ismember(ROIname, format_map.keys)
                format_map(ROIname) = [format_map(ROIname) d_score(f)];
            else
                format_map(ROIname) = d_score(f);
            end
        end

        if ismember(ROIname, mean_ve.keys)
            mean_ve(ROIname) = [mean_ve(ROIname) mean(ve_list)];
        else
            mean_ve(ROIname) = mean(ve_list);
        end
    end
end

% save data
saveFile = fullfile(ExpDir, 'code', 'justin', 'data.mat');
save(saveFile, 'distinctiveness', 'mean_ve');

%% Plotting
% Plot VE across voxel prop for each hemi
y_lh = [];
y_rh = [];
for r = 1:length(ROIs)
    if contains(ROIs{r}, 'lh')
        y_lh = [y_lh mean_ve(ROIs{r})'];
    else
        y_rh = [y_rh mean_ve(ROIs{r})'];
    end
end

savePlot(prop_voxels, y_lh, [0 1], ROI_types, 'Mean VE', 'Mean VE vs. Proportion of Voxels (LH)', ...
    ImageDir, 'lh_ve_plot.png');
savePlot(prop_voxels, y_rh, [0 1], ROI_types, 'Mean VE', 'Mean VE vs. Proportion of Voxels (RH)', ...
    ImageDir, 'rh_ve_plot.png');
fprintf('Saved VE vs. prop vox plots\n');


% Plot distinctiveness across voxel prop for each format
for f = 1:length(formats)
    y_lh = [];
    y_rh = [];
    format = formats{f}
    map = distinctiveness(format);
    for r = 1:length(ROIs)
        if contains(ROIs{r}, 'lh')
            y_lh = [y_lh map(ROIs{r})'];
        else
            y_rh = [y_rh map(ROIs{r})'];
        end
    end
    
    title = ['Distinctiveness (' format ') vs. Proportion of Voxels'];
    filename = [format '_plot.png'];
    savePlot(prop_voxels, y_lh, [-0.1 0.1], ROI_types, 'Distinctiveness', [title ' (LH)'], ...
        ImageDir, ['lh_' filename]);
    savePlot(prop_voxels, y_rh, [-0.1 0.1], ROI_types, 'Distinctiveness', [title ' (RH)'], ...
        ImageDir, ['rh_' filename]);
end
fprintf('Saved distinctiveness vs. prop vox plots\n');


%% Functions
function [] = savePlot(x, y, ylims, leglab, ylab, titlestr, ImageDir, name)
plot(x, y, 'LineWidth', 2);
ylim(ylims);
title(titlestr);
xlabel('Proportion of Voxels');
ylabel(ylab);
legend(leglab, 'Location', 'Northeast');
saveas(gcf, fullfile(ImageDir, name));
close
end

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


% distinctiveness score at category-level
function [score] = distinctiveness_category(rsm)
cat_scores = zeros(1, 8);

for condition=1:8
    start_dyn = ((condition - 1) * 24) + 1;
    start_sta = start_dyn + 192;

    % without cluster
    dyn_to_other = rsm(start_dyn:start_dyn + 23, [1:start_dyn - 1, ...
        start_dyn + 24:start_sta - 1, start_sta:end]);
    sta_to_other = rsm(start_sta:start_sta + 23, [1:start_dyn - 1, ...
        start_dyn + 24:start_sta - 1, start_sta:end]);
    without_cluster = [flat(dyn_to_other) flat(sta_to_other)];
    a = mean(without_cluster);

    % within cluster
    dyn_cluster = aboveDiagonalFlat(rsm(start_dyn:start_dyn + 23, start_dyn:start_dyn + 23));
    sta_cluster = aboveDiagonalFlat(rsm(start_sta:start_sta + 23, start_sta:start_sta + 23));
    across_format_within_cat = rsm(start_dyn:start_dyn + 23, start_sta:start_sta + 23);
    within_cluster = [dyn_cluster sta_cluster flat(across_format_within_cat)];
    b = mean(within_cluster);

    cat_scores(condition) = (b - a)/2;
end

score = mean(cat_scores);

end


% distinctiveness score at animacy-level
function [score] = distinctiveness_animacy(rsm)
ani_to_ani_dyn = aboveDiagonalFlat(rsm(25:120, 25:120));
ani_to_ani_sta = aboveDiagonalFlat(rsm(217:312, 217:312));
ani_to_ani_across = rsm(25:120, 217:312);
within_cluster_ani = [ani_to_ani_dyn ani_to_ani_sta ...
    flat(ani_to_ani_across)];
b_ani = mean(within_cluster_ani);

ina_words = aboveDiagonalFlat(rsm(1:24, 1:24));
ina_mixed = aboveDiagonalFlat(rsm(121:216, 121:216));
ina_sta = aboveDiagonalFlat(rsm(313:end, 313:end));
within_cluster_ina = [ina_words ina_mixed ina_sta ...
    flat(rsm(1:24, [121:216, 313:end])) flat(rsm(121:216, 313:end))];
b_ina = mean(within_cluster_ina);

without_cluster = [flat(rsm(1:24, [25:120, 217:312])) ...
    flat(rsm(25:120, [121:216, 313:end])) ...
    flat(rsm(121:216, 217:312)) ...
    flat(rsm(217:312, 313:end))];
a = mean(without_cluster);

distinctiveness_ani = (b_ani - a)/2;
distinctiveness_ina = (b_ina - a)/2;
score = mean([distinctiveness_ani distinctiveness_ina]);

end


% distinctiveness score at motion-level
function [score] = distinctiveness_motion(rsm)
dyn_to_dyn = mean(aboveDiagonalFlat(rsm(1:192, 1:192)));
sta_to_sta = mean(aboveDiagonalFlat(rsm(193:end, 193:end)));
dyn_to_sta = mean(mean(rsm(1:192, 193:end)));

distinctiveness_dyn = (dyn_to_dyn - dyn_to_sta)/2;
distinctiveness_sta = (sta_to_sta - dyn_to_sta)/2;

score = mean([distinctiveness_dyn distinctiveness_sta]);

end
