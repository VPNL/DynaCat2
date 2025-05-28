fprintf('Initializing Beth Matlab....\n');
% fMRI/dMRI
%------------------ vistasoft ------------------------%
addpath(genpath('/share/kalanit/users/brispoli/matlab/vistasoft'));
%----------------------------------------------------%
%------------------ SPM8 ------------------------%
addpath(genpath('/share/kalanit/software/spm8'));
%-----------------------------------------------------%
%--------------------- dataTools --------------------------------%
dataToolsPath = '/share/kalanit/biac2/kgs/dataTools';
addpath(dataToolsPath);
fprintf('RAID added to work path\n');
%-----------------------------------------------------%
%---------------- VPNL tools ----------------------%
addpath(genpath('/share/kalanit/biac2/kgs/projects/GitHub/VPNLtools'));
%----------------------------------------------------%
%---------------- KNK tools ----------------------%
addpath(genpath('/share/kalanit/software/knkutils'));
%----------------------------------------------------%
% dMRI
%---------------- mrTrix ----------------------%
trixPath = '/usr/lib/mrtrix/bin';
addpath(genpath(trixPath))
fprintf('mrTrix added to work path\n');
%----------------------------------------------------%
%---------------- LiFE/AFQ----------------------%
addpath(genpath('/share/kalanit/software/afq'));
addpath(genpath('/share/kalanit/software/encode'));
%----------------------------------------------------%
% qMRI
%---------------- mrQ  ----------------------%
addpath(genpath('/share/kalanit/software/mrQ'));
%----------------------------------------------------%
% cortical surface tools
%------------ FreeSurfer -----------------------------%
fshome = getenv('FREESURFER_HOME');
fsmatlab = sprintf('%s/matlab',fshome);
if (exist(fsmatlab) == 7)
    addpath(genpath(fsmatlab));
end
clear fshome fsmatlab;
%-----------------------------------------------------%
%------------ FreeSurfer FAST ------------------------%
fsfasthome = getenv('FSFAST_HOME');
fsfasttoolbox = sprintf('%s/toolbox',fsfasthome);
if (exist(fsfasttoolbox) == 7)
    path(path,fsfasttoolbox);
end
clear fsfasthome fsfasttoolbox;
%-----------------------------------------------------%
%-----------------------------------------------------%
%my code
addpath(genpath('/share/kalanit/biac2/kgs/projects/DynaCat/code'));
% addpath(genpath('/share/kalanit/biac2/kgs/projects/DynamicCategories/code/funcitons/'));
%christina's code
addpath(genpath('/share/kalanit/biac2/kgs/projects/bbfloc/code/bbfLoc'));

fprintf('\nAll done!\n');