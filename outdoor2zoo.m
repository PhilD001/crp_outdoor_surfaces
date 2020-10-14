function outdoor2zoo(fld)
% OUTDOOR2ZOO is a custom function to convert data from outdoor data set to zoo format
%
% Arguments
%   fld ... str. Full path to data folder

if nargin==0
    fld=uigetfolder;
end
tic
cd(fld)

% some hardcoded values
side = 'R'; % only extract right side gait cycles
chns = {'Knee', 'Hip'};
sample_rate = 100;

% get files and batch process (each mat file is a separate subject)
fl = engine('fld', fld, 'extension', '.mat');
r = load(fl{1}, '-mat');
r = r.FinalOutput;

cons = fieldnames(r);

for c = 1:length(cons)
    
    if isempty(r.(cons{c}))
        disp(['no data for condition ', cons{c}])
        continue
    end
    
    subs = fieldnames(r.(cons{c}));
    for s = 1:length(subs)
        trials = fieldnames(r.(cons{c}).(subs{s}));
        for t = 1:length(trials)
                   
            % create a zoo file each with 3 gait cycles inside a single "cycle" here
            cycles = fieldnames(r.(cons{c}).(subs{s}).(trials{t}).(side).Knee);
            for n= 1:length(cycles)
                fname = [subs{s}, '_', cons{c}, '_', num2str(trials{t}),'_',cycles{n}, '.zoo'];
                disp(['creating zoo file for ', fname])
                evts =  r.(cons{c}).(subs{s}).(trials{t}).(side).gait.(cycles{n});
                
                data = struct;
                data.zoosystem = setZoosystem(fname);
                data.zoosystem.Units = struct;
                data.zoosystem.Video.Freq = sample_rate;
                data.zoosystem.AVR = 0;
                for ch = 1:length(chns)
                    ndata =  r.(cons{c}).(subs{s}).(trials{t}).(side).(chns{ch}).(cycles{n});
                    data = addchannel_data(data, [side, chns{ch}], ndata, 'video');
                    if ch == 1
                        data.([side, chns{ch}]).event.FS = [evts(1) 0 0];
                        data.([side, chns{ch}]).event.FO = [evts(2) 0 0];
                    end
                end
                
               
                
                % Save all into to file
                zsave(fname,data)        
            end
        end 
    end 
end



%---SHOW END OF PROGRAM-------------------------------------------------------------------------
%
disp(' ')
disp('**********************************')
disp('Finished converting data in: ')
toc
disp('**********************************')

