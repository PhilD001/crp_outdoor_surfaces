%-----------------------------------------------------

% Calibration_code.m (Vaibhav Shah dd/mm/yy)

%--------------------------------------------------------------------------

%


% Hip and Knee flexion-extension angle estimation for angle correction


%

% INPUT

% 1. Filtered data of Acceleration, angular rate and magnetism collected using five IMU sensors.Open-source database by Lou et al. (2020)

% Processing  

% 1. Estimating joint angles. 
% 


% OUTPUT

% 1. Struct named "Calibration" that will have participants wise angle esstimation that will be used for angle correction, and will be saved in the folder as "Calibration_data". 

%--------------------------------------------------------------------------


clc; clear all;
load('data.mat')
warning off
s=1; %ID of the participant
e=30; %Number of Subject/run loops for x times
Counthip =0; % Counter for the output that comes under set condi (error data removed)
Countknee = 0;
Totalknee =0;
Totalhip =0;% Total possible
irowknee  = 1; % Row in the out put Struct
irowhip = 1;% Row in the out put Struct
jcol  =1; % Column in the out put Struct
STable= table();
rowR=1;
rowL=1;
KneeFlx=table();
for s= 1:e  % Loop for the participants
    Sensor = "ID" + s; %ID of dataset
%     s= s+1; % For the next participant data set
    T= 1; % Strat Trial
    Te = 3; % Number of Trial loop you want to run
    Sensor;

    
    for T= 1:Te %Trial loop
        SideSel= 1; %Side select 1= Right and 2 =Left
        Side = ["R";"L"];
        plotColor = ['b';'r']; % Blue for Right; Red for Left
        for SideSel = 1:2 % 1st loop for Right and 2nd for left
            Heelstrike=[];
            Heel=[];
            
            Toe=[];
            % For Data Selection
            shank = "shank"+ Side(SideSel,1);
            thigh = "thigh"+ Side(SideSel,1);
            trunk = "trunk";
            Gmissing= [data.(Sensor).(shank).Gyr_X{T,1} data.(Sensor).(shank).Gyr_Y{T,1} data.(Sensor).(shank).Gyr_Z{T,1}]; % Shank Gyro data
            [Gmissingrow, Gcol] = size(Gmissing);
            Start = 1; %Start of Data row
            m1 = data.(Sensor).(shank).MissingCount{T,1};
            m2 = data.(Sensor).(thigh).MissingCount{T,1};
            m3 = data.(Sensor).(trunk).MissingCount{T,1};
            
            if m1 > m2
                if m3 > m1
                    missing= m3;
                else
                    missing= m1;
                end
            else
                if m3 > m2
                    missing= m3;
                else
                    missing= m2;
                end
            end
            %missing data points
            End = Gmissingrow- missing; %End of Data row
           
               
                
                %Loading data
                Ashank= [data.(Sensor).(shank).Acc_X{T,1}(Start:End) data.(Sensor).(shank).Acc_Y{T,1}(Start:End) data.(Sensor).(shank).Acc_Z{T,1}(Start:End)];
                Gshank= [data.(Sensor).(shank).Gyr_X{T,1}(Start:End) data.(Sensor).(shank).Gyr_Y{T,1}(Start:End) data.(Sensor).(shank).Gyr_Z{T,1}(Start:End)];
                Mshank= [data.(Sensor).(shank).Mag_X{T,1}(Start:End) data.(Sensor).(shank).Mag_Y{T,1}(Start:End) data.(Sensor).(shank).Mag_Z{T,1}(Start:End)];
                Athigh= [data.(Sensor).(thigh).Acc_X{T,1}(Start:End) data.(Sensor).(thigh).Acc_Y{T,1}(Start:End) data.(Sensor).(thigh).Acc_Z{T,1}(Start:End)];
                Gthigh= [data.(Sensor).(thigh).Gyr_X{T,1}(Start:End) data.(Sensor).(thigh).Gyr_Y{T,1}(Start:End) data.(Sensor).(thigh).Gyr_Z{T,1}(Start:End)];
                Mthigh= [data.(Sensor).(thigh).Mag_X{T,1}(Start:End) data.(Sensor).(thigh).Mag_Y{T,1}(Start:End) data.(Sensor).(thigh).Mag_Z{T,1}(Start:End)];
              

                
                %Loading Complete
                % Calling Complementry Filter
                GainA = 0.01;
                GainM = 0.01;
                %                       FUSE=ahrsfilter;
                FUSE = complementaryFilter('AccelerometerGain', GainA,'MagnetometerGain',GainM);
                %Calling complete
                
                % Fusing Acc,Gyro and Mag data with complementry
                % filter and output Quaternion and Anguler velocity
                [Qshank,Angushank] = FUSE(Ashank,Gshank,Mshank);
                %                       FUSE=ahrsfilter;
                FUSE = complementaryFilter('AccelerometerGain', GainA,'MagnetometerGain',GainM);
                
                [Qthigh,Anguthigh] = FUSE(Athigh,Gthigh,Mthigh);
                %if not called two times output will come nha
                

                %Knee flex/ext Angle Start
                
                Qknee =conj(Qthigh) .*Qshank;
                
                Qknee=quat2eul(Qknee);
                
                Knee = -rad2deg(Qknee(:,2));
                
                %Knee flex/ext Angle End

                if SideSel==1
                KneeFlx{rowR,1} = Knee(1,1);
                KneeFlx.Properties.VariableNames{1} = 'Right';

                rowR=rowR+1;
                else
                KneeFlx{rowL,2} = Knee(1,1);  
                KneeFlx.Properties.VariableNames{2} = 'Left';

                rowL=rowL+1;
                end
                
      end
        
        %If condistion for Trials
  
    end
          KneeFlx=table2struct(KneeFlx);
  
            STable(1,s)= table(struct('Knee',KneeFlx));
    STable.Properties.VariableNames{s} = str2mat(Sensor);
    KneeFlx=table();

    rowR=1;
    rowL=1;
end
   Calibration= table2struct(STable);
save('Calibration_data.mat','Calibration')