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

%% STEP 1: CRP process

limb = {'Ipsi', 'Contra'};
group = {'CP', 'TD'};
cond = {'SpinTurn','Straight'};
for g = 1:length(group)
    
    for s = 1:length(cond)
        
        participant = GetSubDirsFirstLevelOnly([fld, filesep, group{g}]);
        for p = 1:length(participant)
            
            fld_gsp = [fld, filesep, group{g}, filesep, participant{p}, filesep, cond{s}];
            fl = engine('fld', fld_gsp, 'extension', 'zoo');
            
            disp(' ')
            disp(['extracting data from ', num2str(length(fl)), ' trials for subject ', ...
                participant{p}, ' on surface ', cond{s}, ' ...'])
            
            % initialize stk for each CRP metric for a given group/surface/participant -------------
            IpsiKHStanceCRP_stk = ones(length(fl), 101);
            ContraKHStanceCRP_stk = ones(length(fl), 101);
            
            IpsiAKStanceCRP_stk = ones(length(fl), 101);
            ContraAKStanceCRP_stk = ones(length(fl), 101);
            
            % loop through all trials for a given group/surface/participant ------------------------
            for f = 1:length(fl)
                [file_pth, file_name, ext] = fileparts(fl{f});
                disp(file_name)
                data = load(fl{f}, '-mat');
                data = data.data;
                
                % Extracts Stance phase indices based on existing gait events ------------
                FSapex = data.SACR_x.event.FSapex(1);      % middle portion
                FSminus1 = data.SACR_x.event.FSminus1(1);  % step before apex
                FSplus1 = data.SACR_x.event.FSplus1(1);    % step after apex
                FSplus2 = data.SACR_x.event.FSplus2(1);    % 2 step after apex
                
                Stance = (FSapex:FSplus1);
                PadStance = (FSminus1:FSplus2); % pads stance phase to next nearest gait event.
                
                for l = 1:length(limb)
                    
                    % Extract Joint Angles ---------------------------------------------------------
                    Hip = data.([limb{l}, 'HipAngles_x']).line;
                    Knee = data.([limb{l}, 'KneeAngles_x']).line;
                    Ankle = data.([limb{l}, 'AnkleAngles_x']).line;
                    
                    % Hip, knee, and ankle phase angle padded to nearest event on either side ------
                    HipCyclePhase  = phase_angle(Hip(PadStance));
                    HipStancePhase = HipCyclePhase(Stance);
                    
                    KneeCyclePhase  = phase_angle(Knee(PadStance));
                    KneeStancePhase = KneeCyclePhase(Stance);
                    
                    AnkleCyclePhase  = phase_angle(Ankle(PadStance));
                    AnkleStancePhase = AnkleCyclePhase(Stance);
                    
                    % CRP calculations -------------------------------------------------------------
                    KHStanceCRP = CRP(KneeStancePhase,HipStancePhase);
                    
                    AKStanceCRP = CRP(AnkleStancePhase, KneeStancePhase);
                    
                    % Time Normalizes CRP curves to 100 percent (101 points) -----------------------
                    KHStanceCRP_Norm = normalize_line(KHStanceCRP, 100, 'spline');
                    AKStanceCRP_Norm = normalize_line(AKStanceCRP, 100, 'spline');
                    
                    if strcmp(limb{l}, 'Ipsi')
                        IpsiKHStanceCRP_stk(f, :) = KHStanceCRP_Norm;
                        IpsiAKStanceCRP_stk(f, :) = AKStanceCRP_Norm;
                        
                    elseif strcmp(limb{l}, 'Contra')
                        ContraKHStanceCRP_stk(f, :) = KHStanceCRP_Norm;
                        ContraAKStanceCRP_stk(f, :) = AKStanceCRP_Norm;
                    end
                end
            end
            
            % delete temp files --------------------------------------------------------------------
            for f = 1:length(fl)
                java.io.File(fl{f}).delete();
            end
            [rws, cols] = size(IpsiKHStanceCRP_stk);
            if rws < 2
                disp('insufficient trials for CRP analysis, setting trial as outlier with 999 code')
                
                IpsiKHStance_MARP = 999*ones(1,101);
                IpsiAKStance_MARP = 999*ones(1,101);
                ContraKHStance_MARP = 999*ones(1,101);
                ContraAKStance_MARP = 999*ones(1,101);
                IpsiKHStance_DP = 999*ones(1,101);
                IpsiAKStance_DP = 999*ones(1,101);
                ContraKHStance_DP = 999*ones(1,101);
                ContraAKStance_DP = 999*ones(1,101);
            else
                % compute mean absolute relative phase (MARP) ------------------------------------------
                IpsiKHStance_MARP = mean(IpsiKHStanceCRP_stk);
                IpsiAKStance_MARP = mean(IpsiAKStanceCRP_stk);
                
                ContraKHStance_MARP = mean(ContraKHStanceCRP_stk);
                ContraAKStance_MARP = mean(ContraAKStanceCRP_stk);
                
                % Compute deviation phase (DP) ---------------------------------------------------------
                IpsiKHStance_DP = std(IpsiKHStanceCRP_stk);
                IpsiAKStance_DP = std(IpsiAKStanceCRP_stk);
                
                ContraKHStance_DP = std(ContraKHStanceCRP_stk);
                ContraAKStance_DP = std(ContraAKStanceCRP_stk);
            end
            
            % create new zoo file for ensembled (mean, std) data -----------------------------------
            fl_ens = [file_pth, filesep, file_name(1:end-2), 'ens', ext];
            data_ens = struct;
            
            % save time series info to file --------------------------------------------------------
            data_ens.IpsiKHStance_MARP.line = IpsiKHStance_MARP;
            data_ens.IpsiAKStance_MARP.line = IpsiAKStance_MARP;
            
            data_ens.ContraKHStance_MARP.line = ContraKHStance_MARP;
            data_ens.ContraAKStance_MARP.line = ContraAKStance_MARP;
            
            data_ens.IpsiKHStance_DP.line = IpsiKHStance_DP;
            data_ens.IpsiAKStance_DP.line = IpsiAKStance_DP;
            
            data_ens.ContraKHStance_DP.line = ContraKHStance_DP;
            data_ens.ContraAKStance_DP.line = ContraAKStance_DP;
            
            % compute events at each phase of gait cycle for each MARP, DP curve -------------------
            ch = fieldnames(data_ens);
            for c = 1:length(ch)
                r = data_ens.(ch{c}).line;
                events = struct;
                
                if strfind(ch{c}, 'Stance')
                    events.IC  = [1,  mean(r(1:4)),    0];
                    events.LR  = [5,  mean(r(5:20)),   0];
                    events.MS  = [21, mean(r(21:50)),  0];
                    events.TS  = [51, mean(r(51:81)),  0];
                    events.PSw = [82, mean(r(82:101)), 0];
                end
                
                data_ens.(ch{c}).event = events;
            end
            
            data_ens.zoosystem = setZoosystem(fl_ens);
            zsave(fl_ens, data_ens);
            
        end
    end
end

%% STEP 2: EXTRACT TRIAL BY TRIAL EVENTS TO SPREADSHEET ====================================================

[~, subjects] = extract_filestruct(fld);
cons = cell(1,4);
count = 1;
for i = 1:2
    for j =1:2
        cons{1, count} = [group{i},filesep, cond{j}];
        count= count +1;
    end
end
lcl_evts = {'IC', 'LR', 'MS', 'TS', 'PSw'};
chns = {'IpsiKHStance_MARP','IpsiAKStance_MARP', 'ContraKHStance_MARP','ContraAKStance_MARP', ...
    'IpsiKHStance_DP','IpsiAKStance_DP', 'ContraKHStance_DP','ContraAKStance_DP'};

eventval('fld', fld, 'dim1', cons, 'dim2', subjects, 'ch', chns, ...
    'localevents', lcl_evts, 'globalevents', 'none');

