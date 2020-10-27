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
data_file = [fld_root, filesep, 'data.zip'];

fld = [fld_root, filesep, 'data'];               % setting path for processed data
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

%% STEP 1: Delete unwanted files ==========================================

word = ["Stair","Slope"];       % to identify all files containing Stair and Slope
fl = engine('fld',fld);

for f = 1:length(fl)
    [file_pth, file_name, ext] = fileparts(fl{f});
        if contains(file_name, word)
           delete([file_name, ext])
           disp(['Deleting ', file_name]) 
        end  
end

%% STEP 2: re-organize data in folders ====================================

surface = {'Grass','FlatEven','CobbleStone'};
fl = engine('fld',fld);

for f = 1:length(fl)
        [file_pth, file_name, ext] = fileparts(fl{f});

    for s = 1:length(surface)
        
        for p = 1:30
            
            if ~exist([fld, filesep, num2str(p), filesep, char(surface(s))],'dir')
                mkdir([fld, filesep, num2str(p), filesep, char(surface(s))])
            end
    
            indx = strfind(file_name, ['ID', num2str(p), '_', char(surface(s))]);
            if indx == 1
                movefile ([file_name, '.zoo'], [fld, filesep, num2str(p), filesep, char(surface(s))])
                disp(file_name)
            end
        end
    end
end
%% STEP 3: CRP Analyses ===================================================

fld_data_proc = [fld_root, filesep, 'processed'];
mkdir(fld_data_proc);
copyfile(fld,fld_data_proc)

surface = {'Grass','FlatEven','CobbleStone'};

for s = 1:length(surface)
    
    participant = GetSubDirsFirstLevelOnly([fld_data_proc, filesep]);
    for p = 1:30
            
        fld_ps = [fld_data_proc, filesep, participant{p}, filesep, surface{s}];
        fl = engine('fld', fld_ps, 'extension', 'zoo');

        disp(' ')
        disp(['extracting data from ', num2str(length(fl)), ' trials for subject ', ...
            participant{p},' ...'])

        % initialize stk for each CRP metric for a given surface/participant -------------
        KHCRP_stk = ones(length(fl), 101);

        % loop through all trials for a given surface/participant ------------------------
        for f = 1:length(fl)
            [file_pth, file_name, ext] = fileparts(fl{f});
            disp(file_name)
            data = load(fl{f}, '-mat');
            data = data.data;

            % Extracts Stance phase indices based on existing gait events ------------
            FS = data.RKnee.event.FS(1);  % Foot Strike
            FO = data.RKnee.event.FO(1);  % Foot Off

            % Extract Joint Angles ---------------------------------------------------------
            Hip = data.RHip.line(FS:FO);
            Knee = data.RKnee.line(FS:FO);

            PaddedHip = data.RHip.line; % Pads with additional cycle on both ends
            PaddedKnee = data.RKnee.line;

            % Hip, knee, and ankle phase angle padded to nearest event on either side ------
            HipCyclePhase  = phase_angle(PaddedHip);
            KneeCyclePhase  = phase_angle(PaddedKnee);

            % CRP calculations -------------------------------------------------------------
            KHCycleCRP = CRP(KneeCyclePhase(FS:FO),HipCyclePhase(FS:FO));

            % Time Normalizes CRP curves to 100 percent (101 points) -----------------------
            KHCycleCRP_Norm = normalize_line(KHCycleCRP, 100, 'spline');
            KHCRP_stk(f,:) = KHCycleCRP_Norm; 
        end

        [rws, cols] = size(KHCRP_stk);
        if rws < 2
            disp('insufficient trials for CRP analysis, setting trial as outlier with 999 code')

            KH_MARP = 999*ones(1,101);
            KH_DP = 999*ones(1,101);

        else
            % compute mean absolute relative phase (MARP) ------------------------------------------
            KH_MARP = mean(KHCRP_stk);

            % Compute deviation phase (DP) ---------------------------------------------------------
            KH_DP = std(KHCRP_stk); 
        end

        % create new zoo file for ensembled (mean, std) data -----------------------------------
        fl_ens = [file_pth, filesep, file_name(1:end-14), 'ens', ext];
        data_ens = struct;

        % save time series info to file --------------------------------------------------------
        data_ens.KH_MARP.line = KH_MARP;
        data_ens.KH_DP.line = KH_DP;

        % compute events at each phase of gait cycle for each MARP, DP curve -------------------
        ch = fieldnames(data_ens);
        for c = 1:length(ch)
            r = data_ens.(ch{c}).line;
            events = struct;
            
                events.IC  = [1,  mean(r(1:2)),    0]; % Stance phase
                events.LR  = [3,  mean(r(3:12)),   0];
                events.MS  = [13, mean(r(13:31)),  0];
                events.TS  = [32, mean(r(32:50)),  0];
                events.PSw = [51, mean(r(51:62)), 0];

                events.ISw = [63,  mean(r(63:75)),    0]; % Swing phase
                events.MSw = [76,  mean(r(76:87)),    0];
                events.TSw = [88,  mean(r(88:100)),    0];

            data_ens.(ch{c}).event = events;
        end

        data_ens.zoosystem = setZoosystem(fl_ens);
        zsave(fl_ens, data_ens);

        % delete temp files --------------------------------------------------------------------
        for f = 1:length(fl)
            java.io.File(fl{f}).delete();
        end

     end
end

%% STEP 4: EXTRACT TRIAL BY TRIAL EVENTS TO SPREADSHEET ====================================================

[~, subjects] = extract_filestruct(fld); % gets subs names

lcl_evts = {'IC', 'LR', 'MS', 'TS', 'PSw','ISw','MSw','TSw'};
chns = {'KH_MARP','KH_DP'};

eventval('fld_data_proc', fld_data_proc, 'dim1', surface, 'dim2', subjects, 'ch', chns, ...
    'localevents', lcl_evts, 'globalevents', 'none');
        
