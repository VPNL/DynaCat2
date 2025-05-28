%% DynaCat + Staticat: Anat ROI RSMs, individual stimuli
%   -Computes multivoxel and plots RSMs using indiviudal stimuli as conditions
%       for each subject and then across subjects within anatomical ROIs.
%   -Choose an ampType. (ampType = 'zscore' or ampType = 'subtractedBetas')
%
% BR, JC

%% Set parameters
subjects = {'AX' 'BR' 'CT' 'DO' 'HC' 'IK' 'KP' 'RU' 'RY' 'VL'};
%ROIs = {'lh_IPS' 'rh_IPS' 'lh_STS' 'rh_STS' 'lh_LOTC' 'rh_LOTC' 'lh_VTC' 'rh_VTC'};
ROIs = {'lh_IPS'};

% set parameters
eventsPerBlock = 4; % length of trial is 4s and the TR is 1s
ampType = 'zscore';
plot_type = 'sym_means';
ve = 15;

ExpDir = '/share/kalanit/biac2/kgs/projects/DynaCat/';
addpath(genpath('/share/kalanit/biac2/kgs/projects/DynaCat/code/scripts/rsms/'))

stimuli_icon_folder = '/share/kalanit/biac2/kgs/projects/DynaCat/stimuli/stimuli_icons/stimuli_icons_lower_res/dynacat_staticat_colored_icons';
stimuli_icon_label = 'exp_colors';
stimuli_icon_size = 0.015;
categeory_color_point_size = 17;

num_all_conds = 384;
parfile_folders = 'individual_stim_parfiles_set1';

category_labels = {'words-d', 'dogs-d', 'bodies-d', 'hands-d', 'faces-d', 'cars-d', 'balls-d', 'scenes-d'...
    'words-s', 'dogs-s', 'bodies-s', 'hands-s', 'faces-s', 'cars-s', 'balls-s', 'scenes-s'};

% load stimuli condition info (for group mds plot labels)
all_stimuli_cond_parfile_path = fullfile(ExpDir, 'code', 'experiment_code', 'DynaCat', 'dynacat_stimuli', 'phase_scrambled_set1', 'all_stimuli_conditions.par');
all_cond_parfile = {};
fid = fopen(all_stimuli_cond_parfile_path, 'r');
data = textscan(fid, '%f %s %s %s %f %f %f', 'Delimiter', '\t');
fclose(fid);
for i=1:length(data{1})
    cond_num = data{1}(i);
    cond = data{2}(i);
    experiment = data{3}(i);
    category = data{4}(i);
    color_1 = data{5}(i);
    color_2 = data{6}(i);
    color_3 = data{7}(i);
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

%% Compute + plot MVP, RSM per subject + across subjects
for r = 1:length(ROIs) % loop over ROIs
    ROIname = ROIs{r}

    % initialize across-subject rsm
    across_subject_rsm = nan(num_all_conds, num_all_conds, length(subjects));
    group_distance_matrix = nan(num_all_conds, num_all_conds, length(subjects));

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
        % filtered_data = sorted_parfile_data(~strcmp(sorted_parfile_data(:, 4), 'baseline') & ~strcmp(sorted_parfile_data(:, 4), 'task-repeat'), :);
        cond_nums = cell2mat(sorted_parfile_data(:, 3));
        rowsToKeep = (cond_nums ~= 0) & (cond_nums ~= 999); % conds not seen by others for BR and IK are coded as 999
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

        labels = filtered_parfile_data(:, 4)';
        label_colors = filtered_parfile_data(:, 5);

        % make an array of point shapes for the mds
        point_shapes = cell(size(labels));
        for i = 1:numel(labels)
            if contains(labels{i}, 'dynacat')
                point_shapes{i} = 'o';
            elseif contains(labels{i}, 'staticat')
                point_shapes{i} = '^';
            end
        end

        % numconds = length(labels);
        parfile_cond_nums = cell2mat(filtered_parfile_data(:, 3));
        numconds = length(unique(parfile_cond_nums));

        % ===== plot rsm =====

        % copy ROI over from mrVistaROI folder if nesecary
        ROI_filepath = fullfile(session_path, '3DAnatomy', 'ROIs', [ROIname '.mat']);
        ROI_filepath_in_mrVistaROI_dir = fullfile(session_path, '3DAnatomy', 'mrVistaROIs', [ROIname '.mat']);
        if ~exist(ROI_filepath, 'file')
            copyfile(ROI_filepath_in_mrVistaROI_dir, ROI_filepath);
        end

        % Set vANATOMY
        hg=initHiddenGray(datatype, scan, ROIname);

        % add individual stim parfiles
        runs = list_runs;
        hg = er_assignParfilesToScans(hg, runs, newscans_parfiles);

        % Initialize multivoxel pattern
        mv=mv_init(hg,ROIname,runs,datatype);
        mvs = {mv};

        % Set parameters
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

        % Apply GLM to multivoxel pattern
        mv = mv_applyGlm(mv);

        % get betas for each run and voxel
        amps = mv_amps(mv);

        % removes task condition from amps
        amps = amps(:, 1:numconds);

        % values to plot
        z_values = amps;

        if strcmp(ampType, 'subtractedBetas') || strcmp(ampType, 'zscore')
            % subtract mean beta from each voxel
            mean_amps=mean(amps,2); % calculate the mean for each voxel
            meanmat=mean_amps*ones(1,numconds);
            subbetas=amps-meanmat;
            % the values used for the correlations
            z_values=subbetas;
        end

        if strcmp(ampType, 'zscore')
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

        % Titrate
        voxels_ve = find(mv.glm.varianceExplained > ve/100);
        %[B, voxels] = maxk(mv.glm.varianceExplained, 1000);
        z_values = z_values(voxels_ve',:);

        % Compute crosscorrelation matrices
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

        % save subject rsm data
        subjectDataDir = fullfile(DataDir, subject);
        if ~exist(subjectDataDir, 'dir')
            mkdir(subjectDataDir);
        end
        rsm_filepath = fullfile(subjectDataDir, [ROIname '_' num2str(ve) '_ve_rsm.mat']);
        save(rsm_filepath, 'rsm');
        % save subject all cond rsm data
        rsm_all_cond_filepath = fullfile(subjectDataDir, [ROIname '_' num2str(ve) '_ve_rsm_all_cond.mat']);
        save(rsm_all_cond_filepath, 'all_cond_rsm');

        across_subject_rsm(:, :, s) = all_cond_rsm; % add this subject to across subject rsm

        close all % closing things to save memory

        % clear everything to save memory
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
    mean_rsm = nanmean(across_subject_rsm, 3); % should be 16 by 16 matirx

    % plot across-subject rsm
    figure('color', [1 1 1], 'units','norm', 'position', [0.1 .1 .6 .6]);
    imagesc(mean_rsm, [-.7 .7]); axis('image');
    cmap = mrvColorMaps('coolhot'); colormap(cmap);
    hcb=colorbar; title(hcb, 'correlation');

    plot_title = [ROIname ' across-subjects rsm (n=' num2str(length(subjects)) ...
        ', ve=' num2str(ve) ')'];

    title(plot_title, 'Interpreter','none', 'FontSize',24);
    set(gca,'Xtick', [1:24:max(num_all_conds)], 'XtickLabel',category_labels,'FontSize',10);
    set(gca,'Ytick', [1:24:max(num_all_conds)], 'YtickLabel',category_labels, 'FontSize',10);
    clim([-0.35, 0.35]);
    brighten(0.6)

    % save png
    RSM_filename_png = [ROIname '_' num2str(ve) '_ve.png'];
    outfile = fullfile(DataDir , RSM_filename_png);
    saveas(gcf, outfile);

    % save mean rsm data
    mean_rsm_filepath = fullfile(subjectDataDir, [ROIname '_' num2str(ve) '_ve_rsm.mat']);
    save(mean_rsm_filepath, 'mean_rsm');

    close all 
    clear mean_rsm
end

fprintf('\nAll done!');