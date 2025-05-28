% DynaCat + Staticat: Anat ROI RSMs
% 
%   -Computes multivoxel and plots RSMs for each subject and then across 
%       subjects for anatomical ROIs.
%   -Can do split-half where each side of the matrix is either odd or even
%       runs and the diagonal is the correlation for the same conditions
%       across runs. (plot_type = 'split-half')
%   -Can do mean rsm where both sides are symmetrical, which is what you
%       want if you are going to make an MDS later on. (plot_type = 'sym_means')
%   -Choose an ampType. (ampType = 'zscore' or ampType = 'subtractedBetas')
% 
% BR + KGS 
% 9 May 2024
% JC
% 12 Oct 2024

%% set parameters

subjects = {'AX', 'BR', 'DO', 'HC', 'IK', 'KP', 'RU', 'RY', 'VL'};
ROI_types = {'mtg_bodies' 'itg_bodies' 'ios_words' 'pots_words' ...
    'mog_scenes', 'ips_scenes'};
ROIs = cell(1, length(ROI_types)*2);
for r = 1:length(ROI_types)
    ROIs{r*2 - 1} = strcat('lh_', ROI_types{r});
    ROIs{r*2} = strcat('rh_', ROI_types{r});
end

% Set parameters
numconds = 16; % not including task
eventsPerBlock=4; % length of trial is 4s and the TR is 1s
nrRuns = 8;
ampType = 'zscore';
plot_type = 'sym_means'; % plot_type = 'split-half', plot_type = 'sym_means'

ExpDir='/share/kalanit/biac2/kgs/projects/DynaCat/';
addpath(genpath('/share/kalanit/biac2/kgs/projects/DynaCat/code/scripts/rsms/'))

labels = {'words-d', 'dogs-d', 'bodies-d', 'hands-d', 'faces-d', 'cars-d', 'balls-d', 'scenes-d'...
    'words-s', 'dogs-s', 'bodies-s', 'hands-s', 'faces-s', 'cars-s', 'balls-s', 'scenes-s'};

%% compute and plot mvp and rsm per subject and across subjects
for r=1:length(ROIs) % loop over ROIs
    ROIname = ROIs{r};
    across_subject_rsm = zeros(numconds, numconds, length(subjects));
    
    for s=1:length(subjects) % loop over sessions
        % get subject session info
        subject = subjects{s};
        [session_path, list_runs, datatype, scan] = dynacat_staticat_sessions(subject);
        cd(session_path)
        
        % init session
        hg=initHiddenGray(datatype, scan, ROIname); % had to copy ROI over from mrVistaROIs to ROIs in 3dAnat
        currentview = 'gray';

        if ~isempty(hg.ROIs) % ROI exists
            % load correct parfiles
            parfiles_folder = fullfile(session_path, 'parfiles');
            parfiles_folder = fullfile(parfiles_folder, '*.par');
            parfiles_dir = dir(parfiles_folder);
            parfiles = fullfile({parfiles_dir.folder}, {parfiles_dir.name});
            hg = er_assignParfilesToScans(hg, list_runs, parfiles);
            
            if strcmp(plot_type, 'sym_means')
                runs = list_runs;
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
                amps=mv_amps(mv);
        
                % removes task condition from amps 
                amps=amps(:, 1:numconds);
        
                % values to plot
                z_values=mv_amps(mv);
                % remove task
                z_values=z_values(:, 1:numconds);
                
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
    
                % plot and save subject MVP
                figure('color', [ 1 1 1], 'name', ['MVP' ROIs{r} '_' subjects{s} '_' ampType], 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.4]);
                
                allvals=[z_values];
                maxvalue=max(max(allvals));
                minvalue=min(min(allvals));
                
                imagesc(z_values,[minvalue maxvalue]); colorbar; %set both graphs to have same min and max values
                set(gca,'Xtick', [1:1:numconds], 'XtickLabel',labels,'FontSize',8);
                ylabel('voxel number'); xlabel('condition');
                myTitle = sprintf('MVP across runs %s %s %s', subject, ROIs{r}, ampType);
                title(myTitle, 'Interpreter','none', 'FontSize',10);
                hcb=colorbar; title(hcb, ampType); % adds title to the colorbar
               
                OutDir = fullfile(ExpDir, 'code', 'scripts', 'rsms', 'justin', subject);
                if ~exist(OutDir,'dir')
                    mkdir(OutDir);
                end
                MV_filename_png=[ROIname '_' subject '_mvp_' plot_type '.png'];
                outfile=fullfile(OutDir , MV_filename_png);
                saveas(gcf, outfile, 'png')
        
                % Compute crosscorrelation matrices
                rsm = zeros(numconds);
                for i=1:numconds % calculate correlation between each condition
                    for j=1:numconds
                        c = corrcoef(z_values(:,i),z_values(:,j));
                        rsm(i, j) = c(1, 2);
                    end
                end
    
                across_subject_rsm(:, :, s) = rsm; % add this subject to across subject rsm
            end
    
            % visualize & save individual RSMs as images
            figure('color', [ 1 1 1], 'name', [subject ' '  ROIname],'units','norm', 'position', [ 0.1 .1 .6 .6]);
            imagesc(rsm, [-.7 .7]); axis('image');
            cmap=mrvColorMaps('coolhot'); colormap(cmap);
            hcb=colorbar; title(hcb, 'correlation');
            set(gca,'Fontsize',12);
            myTitle = sprintf('RSM %s %s %s', subject, ROIs{r}, ampType);
            ylabel('Odd Runs'); xlabel('Even Runs')
            
            title(myTitle, 'Interpreter','none', 'FontSize',14);
            set(gca,'Xtick', [1:1:numconds], 'XtickLabel',labels,'FontSize',12);
            set(gca,'Ytick', [1:1:numconds], 'YtickLabel',labels, 'FontSize',12);
            
            OutDir = fullfile(ExpDir, 'code', 'scripts', 'rsms', 'justin', subject);
            if ~exist(OutDir,'dir')
                mkdir(OutDir);
            end
    
            RSM_filename_png=[ROIname '_' subject '_rsm_' plot_type '.png'];
            outfile=fullfile(OutDir , RSM_filename_png);
            saveas(gcf, outfile, 'png')
    
            % save subject rsm data
            results_folder = fullfile(OutDir, 'results');
            if ~exist(results_folder,'dir')
                mkdir(results_folder);
            end
            rsm_filepath = fullfile(results_folder, [ROIname '_mean_rsm.mat']);
            save(rsm_filepath, 'rsm')
            % save mv data struc
            mv_filepath = fullfile(results_folder, [ROIname '_mv.mat']);
            save(mv_filepath, 'mv')
        end
    end

    % get mean across subjects rsms per roi
    rsm_idx = any(across_subject_rsm ~= 0,[1 2]);
    across_subject_rsm = across_subject_rsm(:,:,rsm_idx); % remove subjects without the roi
    mean_rsm = mean(across_subject_rsm, 3); % should be 16 by 16 matirx

    if size(across_subject_rsm, 3) > 0 % if any have roi
        % plot across-subject rsm
        figure('color', [ 1 1 1], 'name', [ROIname ' across-subjects (n=' num2str(size(across_subject_rsm, 3)) ') ' plot_type ': ' ROIname],'units','norm', 'position', [ 0.1 .1 .6 .6]);
        imagesc(mean_rsm, [-.7 .7]); axis('image');
        cmap=mrvColorMaps('coolhot'); colormap(cmap);
        hcb=colorbar; title(hcb, 'correlation');
        set(gca,'Fontsize',12);
        plot_title = [ROIname ' across-subjects rsm (n=' num2str(length(subjects)) ')'];
        ylabel('Odd Runs'); xlabel('Even Runs')
        
        title(plot_title, 'Interpreter','none', 'FontSize',14);
        set(gca,'Xtick', [1:1:numconds], 'XtickLabel',labels,'FontSize',12);
        set(gca,'Ytick', [1:1:numconds], 'YtickLabel',labels, 'FontSize',12);
        % set x and y tick with labels for conditions
        
        OutDir = fullfile(ExpDir, 'code', 'scripts', 'rsms', 'justin');
        if ~exist(OutDir,'dir')
            mkdir(OutDir);
        end
    
        RSM_filename_png=[ROIname '_rsm_across-subjects_ ' plot_type '.png'];
        outfile=fullfile(OutDir , RSM_filename_png);
        saveas(gcf, outfile, 'png')
    
        % save across subject rsm data
        % mean_rsm_filepath = fullfile(ExpDir, 'results', 'anat_rois', [ROIname '_' plot_type '_mean_rsm.mat']);
        save(OutDir, 'mean_rsm')
    end

    close all

end


%% save ROI grouped plots
% to get hemispheres on the correct side, you should list lh first per ROI 
% in your ROI list
%{
% save grouped plots per subject
for s=1:length(subjects)
    subject = subjects{s};
    [session_path, list_runs, datatype, scan] = dynacat_staticat_sessions(subject);
    cd(fullfile(session_path, 'results', 'anat_rois_rsms'))

    % initalize plot
    figure('color', [ 1 1 1], 'name', ['Anat_ROIs_' subjects{s} '_' plot_type], 'Units', 'normalized', 'Position', [0.1, 0.1, 0.5, 1]);
    num_subplot_rows = length(ROIs)/2;

    for r=1:length(ROIs) % loop over ROIs
        ROI = ROIs{r};
        % load roi plot data
        roi_name = [ROI '_mean_rsm.mat'];
        roi_rsm = load(roi_name);
        roi_rsm = roi_rsm.rsm;

        subplot (num_subplot_rows,2,r);
        imagesc(roi_rsm, [-.7 .7]); axis('image');
        cmap=mrvColorMaps('coolhot'); colormap(cmap);
        hcb=colorbar; title(hcb, 'correlation');
        set(gca,'Fontsize',12);
        plot_title = [subject ' ' ROI ' rsm'];
        title(plot_title, 'Interpreter','none', 'FontSize',10);
        set(gca,'Xtick', [1:1:numconds], 'XtickLabel',labels,'FontSize',8);
        set(gca,'Ytick', [1:1:numconds], 'YtickLabel',labels, 'FontSize',8);

    end
    % save figure
    OutDir = fullfile(ExpDir, 'code', 'scripts', 'rsms', 'justin');
    if ~exist(OutDir,'dir')
        mkdir(OutDir);
    end
    RSM_filename_png=[subject '_rsm_per_roi_' plot_type '.png'];
    outfile=fullfile(OutDir , RSM_filename_png);
    saveas(gcf, outfile, 'png')
    % also save to group results folder
    per_subject_rsm_folder = fullfile(ExpDir, 'results', 'anat_rois', 'per_subject_rsm');
    if ~exist(per_subject_rsm_folder,'dir')
        mkdir(per_subject_rsm_folder);
    end
    outfile=fullfile(per_subject_rsm_folder , RSM_filename_png);
    saveas(gcf, outfile, 'png')

end

%% save across-subject grouped plots
cd(fullfile(ExpDir, 'results', 'anat_rois'))
% initalize plot
figure('color', [ 1 1 1], 'name', ['Anat_ROIs_Across-Subjects_' plot_type], 'Units', 'normalized', 'Position', [0.1, 0.1, 1, 1.5]);
num_subplot_rows = length(ROIs)/2;

for r=1:length(ROIs)
    ROI = ROIs{r};
    % load roi plot data
    roi_name = [ROI '_sym_means_mean_rsm.mat'];
    roi_rsm = load(roi_name);
    roi_rsm = roi_rsm.mean_rsm;

    subplot (num_subplot_rows,2,r);
    imagesc(roi_rsm, [-.7 .7]); axis('image');
    cmap=mrvColorMaps('coolhot'); colormap(cmap);
    hcb=colorbar; title(hcb, 'correlation');
    set(gca,'Fontsize',12);
    plot_title = ['across-subjects ' ROI ' rsm'];
    title(plot_title, 'Interpreter','none', 'FontSize',10);
    set(gca,'Xtick', [1:1:numconds], 'XtickLabel',labels,'FontSize',8);
    set(gca,'Ytick', [1:1:numconds], 'YtickLabel',labels, 'FontSize',8);
end
mainTitle = ['across-subjects rsms (n=' num2str(length(subjects)) ')'];
sgtitle(mainTitle);

%OutDir= fullfile(ExpDir, 'results', 'anat_rois');
if ~exist(OutDir,'dir')
    mkdir(OutDir);
end
RSM_filename_png=sprintf('across_subjects_RSMs_%s.png', plot_type');
outfile=fullfile(OutDir , RSM_filename_png);
saveas(gcf, outfile, 'png')
%}
disp('done')
