% Master script to process data for continuous relative phase analysis of
% outdoor surface gait using IMUs public database 
% (see : https://doi.org/10.1038/s41597-020-0563-y)
%
% Instructions:
% 1) Add biomechZoo and current repo to your matlab path
% 2) Change current working directory to root of crp_outdoor_surfaces directory
% 4) run this script for crp processing
% 5) see the excel sheet generated in /Statistics/eventval.xls for MARP and DP values that
%    can be used for statistical analyses (TO BE CODED)

% Notes
% - tested on biomechZoo v1.4.8 (type zooinfo) under Matlab 2020a

%% STEP 0: prepare data

fld_root = pwd;                                  % make sure pwd points to root of repo
data_file = [fld_root, filesep, 'test_data.zip'];

fld = [fld_root, filesep, 'test_data'];               % setting path for processed data
fld_stats = [fld_root, filesep, 'Statistics'];
if exist(fld, 'dir')
    disp('removing old processed data folder...')
    rmdir(fld, 's')
end
disp(['unzipping data file ', data_file])
f = unzip(data_file);


if exist(fld_stats, 'dir')
    disp('removing old stats folder...')
    rmdir(fld_stats, 's')
end
disp('creating folder for excel sheet output...')
mkdir(fld_stats);        % create empty folder for eventval.xls

outdoor2zoo(fld)

%% STEP 1: CRP analysis



