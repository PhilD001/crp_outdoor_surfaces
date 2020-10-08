%-----------------------------------------------------

% Knee_Hip_Flex_Ext_Estimation.m (Vaibhav Shah dd/mm/yy)

%--------------------------------------------------------------------------

%Run program for 7 times, one time for each surface
% Variable "select" will let you let you choose the surface

% Hip and Knee flexion-extension angle estimation in combination with gait phase detection algorithm


%

% INPUT

% 1. Filtered data of Acceleration, angular rate and magnetism collected using five IMU sensors.Open-source database by Lou et al. (2020)

% Processing  

% 1. Shank sensor’s angular velocity for estimating heel strikes.
% 2. Estimating joint angles using complimentary filter-based sensor fusion algorithm.
% 3. Angle correction


% OUTPUT

% 1. Struct named "FinalOutput" that will have participants wise angle esstimation, which will be saved in the folder as "FinalOutput.mat". 
% Plot of both the joint angles of all the participants in a single 2X2  figure. Blue = right side; Red = Left side.

%--------------------------------------------------------------------------
%%
% clc; clear all;
% 
% load('data.mat');
% load('Calibration_data.mat');
%%  Uncomment only when running for the first time
% FinalOutput= struct('FlatEven',[],'CobbleStone',[],'StairUp',[],'StairDown',[],'SlopeUp',[],'SlopeDown',[],'Grass',[]);
%% uncomment after running the program 1 time
%load('FinalOutput.mat')
%% Variables that should be set
select=1; % 1=FlatEven;2=CobbleStone;3=StairUp;4=StairDown;5=SlopeUp;6=SlopeDown;7=Grass;
s=1; %ID of the participant
e=30; %Number of Subject/run loops for x times; e=30 will run for all 30 participants 
Te = 6; % Number of trial loop you want to run--Shoul be set 6 as there are six numbers of trial on each surfaces
Sideselect= 2; % For both side=2

%% Counter for the output that comes under set condition  (error data removed)
Counthip =0;
Count = 0;
Total =0;
Totalhip =0;
% Total possible
outrow  = 1; % Row in the output Struct of each trial
irowhip = 1;% Row in the output Struct of each trial
STable= table(); % participant output table
gait_cycle= table();
    KneeFlx = table();
    HipFlx = table();
SideTable= table();
SurfaceSelect= ["FlatEven","CobbleStone","StairUp","StairDown","SlopeUp","SlopeDown","Grass"];
Tselect=[4,10,16,17,28,29,52];

%%
for s= 1:e  % Loop for the participants
    %% participant select
    Sensor = "ID" + s; %ID of dataset
    %% Selecting calibration data
    A = table2array(struct2table(Calibration.(Sensor).Knee));
    %% Trial Select
    T= Tselect(select); % Strat trial number according to surface
    TrialNumber=1;
    TrialTable= table();
    %set peak limit according to surface
    %--------------------------------------------------------------------------------------------------
    %Surface |FlatEven|ColbbeStone|StairUP          |     StairDown   |SlopeUp          |SlopeDown        |Grass|%%
    %PeakLim |0.7     |  0.7      |        0.52     |     0.6         |   0.67          |    0.6          | 0.7 |%%
    %  T     |4-9     |  10-15    |16-18-20-22-24-26|17-19-21-23-25-27|28-30-32-34-36-38|29-31-33-35-37-39|52-57|%%
    %-------------------------------------------------------------------------------------------------
    %% peak detection limit select
    if T <16
        PeakLim=0.7; %% FlatEven and ColbbeStone
    elseif T<28
        if T > 15
            even = rem(T,2);
            if even ==0
                if s==29
                    PeakLim=0.4; %StairUP(only for participant 29)
                else
                    PeakLim=0.52;%StairUP
                end
            else
                PeakLim=0.6; %StairDown
            end
        end
    elseif T> 27
        if T < 40
            even = rem(T,2);
            if even ==0
                PeakLim=0.67; % %SlopeUp
            else
                PeakLim=0.6; %SlopeDown
            end
        end
        
    elseif T>51
        PeakLim=0.7; %Grass
    end
    %% Loop for trials
    for Trial= 1:Te %Trial loop (Trial can not be set as T because for some we need to take T+1 and for some we need T+2 trail number)
        
        
        Side = ["R";"L"]; %Side select
        plotColor = ['b';'r']; % Blue for Right; Red for Left
        %% Loop for side
        for SideSel = 1: Sideselect %Side select 1= Right and 2 =Left
            Heelstrike=[];
            Heel=[];
            GyroPeaks=[];
            %% Data Selection
            shank = "shank"+ Side(SideSel,1);
            thigh = "thigh"+ Side(SideSel,1);
            trunk = "trunk";
            %% End number of Row selection
            Gmissing= [data.(Sensor).(shank).Gyr_X{T,1} data.(Sensor).(shank).Gyr_Y{T,1} data.(Sensor).(shank).Gyr_Z{T,1}]; % Shank Gyro data
            [Gmissingrow, Gcol] = size(Gmissing);
            Start = 1; %Start of Data row
            %missing value conunt select form the database
            m1 = data.(Sensor).(shank).MissingCount{T,1};
            m2 = data.(Sensor).(thigh).MissingCount{T,1};
            m3 = data.(Sensor).(trunk).MissingCount{T,1};
            
            if m1 > m2
                if m3 > m1
                    missing= m3;
                else
                    missing= m1;
                end
            elseif m3 > m2
                missing= m3;
            else
                missing= m2;
            end
            
            %missing data points
            End = Gmissingrow- missing; %End row
            
            if End >100 % Only if the total rows are higher than 100
                %% Start of Heel strike detection
                Gs= [data.(Sensor).(shank).Gyr_X{T,1}(Start:End) data.(Sensor).(shank).Gyr_Y{T,1}(Start:End) data.(Sensor).(shank).Gyr_Z{T,1}(Start:End)]; % Shank Gyro data
                G=sqrt(Gs(:,1).^2+Gs(:,2).^2+Gs(:,3).^2); %Normalized Gyro Data
                
                [GRow, Gcol] = size(G); % Count Row and Columm of Gyro Shank
                limitPeakH = PeakLim * max(G,[],'all'); %Set limit peak (X to Max value of Gyro)
                limitPeakD = 50; %Min Distance between two peaks
                % Peak ditection Start
                [~,GyroPeaks,~,~] = findpeaks(G,[1:GRow],...
                    'MinPeakProminence',limitPeakH,...
                    'MinPeakDistance',limitPeakD);
                
                GyroPeaks= GyroPeaks.';
                [iend, row] = size(GyroPeaks);
                
                for i= 1:iend-1
                    [TF,~]=islocalmin(G(GyroPeaks(i):GyroPeaks(i+1)));
                    TFmin=find(TF==1);
                    Heel(i)= GyroPeaks(i)+TFmin(1)-1;
                end
                Heelstrike=Heel;
                [Heelstrikec, Heelstrikei] = size(Heelstrike); %Count number of Row and Column
                
                %% Knee and Hip angle estimation
                
                %Loading data
                Ashank= [data.(Sensor).(shank).Acc_X{T,1}(Start:End) data.(Sensor).(shank).Acc_Y{T,1}(Start:End) data.(Sensor).(shank).Acc_Z{T,1}(Start:End)];
                Gshank= [data.(Sensor).(shank).Gyr_X{T,1}(Start:End) data.(Sensor).(shank).Gyr_Y{T,1}(Start:End) data.(Sensor).(shank).Gyr_Z{T,1}(Start:End)];
                Mshank= [data.(Sensor).(shank).Mag_X{T,1}(Start:End) data.(Sensor).(shank).Mag_Y{T,1}(Start:End) data.(Sensor).(shank).Mag_Z{T,1}(Start:End)];
                
                Athigh= [data.(Sensor).(thigh).Acc_X{T,1}(Start:End) data.(Sensor).(thigh).Acc_Y{T,1}(Start:End) data.(Sensor).(thigh).Acc_Z{T,1}(Start:End)];
                Gthigh= [data.(Sensor).(thigh).Gyr_X{T,1}(Start:End) data.(Sensor).(thigh).Gyr_Y{T,1}(Start:End) data.(Sensor).(thigh).Gyr_Z{T,1}(Start:End)];
                Mthigh= [data.(Sensor).(thigh).Mag_X{T,1}(Start:End) data.(Sensor).(thigh).Mag_Y{T,1}(Start:End) data.(Sensor).(thigh).Mag_Z{T,1}(Start:End)];
                
                Atrunk= [data.(Sensor).(trunk).Acc_X{T,1}(Start:End) data.(Sensor).(trunk).Acc_Y{T,1}(Start:End) data.(Sensor).(trunk).Acc_Z{T,1}(Start:End)];
                Gtrunk= [data.(Sensor).(trunk).Gyr_X{T,1}(Start:End) data.(Sensor).(trunk).Gyr_Y{T,1}(Start:End) data.(Sensor).(trunk).Gyr_Z{T,1}(Start:End)];
                Mtrunk= [data.(Sensor).(trunk).Mag_X{T,1}(Start:End) data.(Sensor).(trunk).Mag_Y{T,1}(Start:End) data.(Sensor).(trunk).Mag_Z{T,1}(Start:End)];
                
                %Data loaded
                
                %% Fusing Acc,Gyro and Mag data with complementry filter
                
                %  Complementry Filter
                GainA = 0.01;
                GainM = 0.01;
                FUSE = complementaryFilter('AccelerometerGain', GainA,'MagnetometerGain',GainM);
                [Qshank,Angushank] = FUSE(Ashank,Gshank,Mshank);
                
                FUSE = complementaryFilter('AccelerometerGain', GainA,'MagnetometerGain',GainM);
                [Qthigh,Anguthigh] = FUSE(Athigh,Gthigh,Mthigh);
                %if not called three times output will co
                FUSE = complementaryFilter('AccelerometerGain', GainA,'MagnetometerGain',GainM);
                [Qtrunk1,Angtrunk] = FUSE(Atrunk,Gtrunk,Mtrunk);
                %% Trunk rotation
                quat = quaternion([0,0,180],'eulerd','XYZ','point');
                Q= compact(Qtrunk1);
                v=(Q(:,2:4));
                Scal=(Q(:,1));
                V=rotateframe(quat, v);
                V2=[Scal,V];
                Qtrunk = quaternion(V2);
                %% Knee flex/ext Angle Start
                
                Qknee =conj(Qthigh) .*Qshank;
                Qknee=quat2eul(Qknee);
                Calib= min(A(:,SideSel));
                Knee = (-rad2deg(Qknee(:,2)))-Calib;
                
                %Knee flex/ext Angle End
                
                %% Hip flex/ext Angle Start
                Qhip = conj(Qtrunk) .* Qthigh;
                Qhip=quat2eul(Qhip);
                Hip = rad2deg(Qhip(:,2));
                %Hip flex/ext Angle End
                %% Output
                
                for i= 1:Heelstrikei-3
                    OutputX= Heelstrike(i); %ouput Row Start
                    OutputY= Heelstrike(i+3); %output Row end
                    Kneeflx = Knee(OutputX:OutputY);
                    Hipflx  = Hip (OutputX:OutputY);
                    color = plotColor(SideSel,1);
                    Opphip = Hipflx(1,1);
                    Knee_start = Kneeflx(1,1); 
                    if select==1
                                if Opphip <-10
                                    Hipflx= -Hipflx;
                                end
                    end

                    KC = Heelstrike(i+2)-Heelstrike(i+1);
                    
                    if KC < 150
                        subplot(2,1,1)
                        plot (Kneeflx, color);
                        hold on
                        subplot(2,1,2)
                        plot (Hipflx, color);
                        hold on
                        Count= Count+1;
                        if outrow==1

                            KneeFlx(1,1) = {Kneeflx};
                            HipFlx(1,1) = {Hipflx};
                            gait_cycle(1,outrow) = {[(Heelstrike(i+1)-Heelstrike(i));(((Heelstrike(i+1)-Heelstrike(i))+(Heelstrike(i+2)-Heelstrike(i+1))))]};
                            cycle= "cycle_" + (outrow);
                            gait_cycle.Properties.VariableNames{outrow}= str2mat(cycle);
                            KneeFlx.Properties.VariableNames{outrow}= str2mat(cycle);
                            HipFlx.Properties.VariableNames{outrow}= str2mat(cycle);
                            
                           outrow= outrow+1;
                            
                        else
                            KneeFlx(1,outrow) = {Kneeflx};
                            HipFlx(1,outrow) = {Hipflx};
                            gait_cycle(1,outrow) = {[(Heelstrike(i+1)-Heelstrike(i));(((Heelstrike(i+1)-Heelstrike(i))+(Heelstrike(i+2)-Heelstrike(i+1))))]};
                            cycle= "cycle_" + (outrow);
                            gait_cycle.Properties.VariableNames{outrow}= str2mat(cycle);
                            KneeFlx.Properties.VariableNames{outrow}= str2mat(cycle);
                            HipFlx.Properties.VariableNames{outrow}= str2mat(cycle);
                            
                            outrow= outrow+1;
                        end
                        
                    else
                       
                    end
                    Total = Total+1;

                end

            end
    SideTable(1,SideSel)= table(struct('Knee',table2struct(KneeFlx),'Hip',table2struct(HipFlx),'gait',table2struct(gait_cycle)));
    SideTable.Properties.VariableNames{SideSel}= str2mat(Side(SideSel,1));
    KneeFlx = table();
    HipFlx = table();
    gait_cycle= table();
    outrow=1;
        end
      
       %% %Next trial selection
        
        if T < 16
            T= T+1;

        elseif T<52
                T= T+2;
            
        else
                T= T+1;
        end
    TrialTable(1,TrialNumber)= table(table2struct(SideTable));
    Pro_trial= "trial" + (TrialNumber);
    TrialTable.Properties.VariableNames{TrialNumber}= str2mat(Pro_trial);
    TrialNumber=TrialNumber+1;
  
    end

   
       
    STable(1,s)= table(table2struct(TrialTable));
    STable.Properties.VariableNames{s} = str2mat(Sensor);
    TrialTable= table();
  
    hold on
end

SucRate = 100*Count/Total;

% 
FinalOutput.(SurfaceSelect(select))=table2struct(STable);
save('FinalOutput.mat','FinalOutput')