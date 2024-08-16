% This parameter file is used for BEHAVIOR task test on Hanlei's PC

warning off;
SetupRand;
HideCursor;
%set_test_gamma;

Param = struct;


%*************************************
%*************Open Window*************
%************************************* 
screens = Screen('Screens');
screenNumber = max(screens);	

white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
gray = round((white+black)/2);
Screen('Preference', 'SkipSyncTests', 0);
[wnd, wRect] = Screen('OpenWindow',screenNumber,gray); 
Screen('BlendFunction', wnd, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

FrameRate = Screen('FrameRate',wnd);
NominalFrameRate = Screen('NominalFrameRate',wnd); 
RefreshDur = Screen('GetFlipInterval',wnd);
Slack = RefreshDur / 2;

%*************************************
%***********Basic Settings************
%*************************************
Param.Settings.ViewDistance   = 400;            % mm 
Param.Settings.ScrnResolution = wRect;
Param.Settings.SquareSize     = wRect(3);       % pixels
Param.Settings.SquareLength   = 700;            % mm 
Param.Settings.PixelPerDegree = 2*Param.Settings.ViewDistance*tan(1/2/180*pi)*Param.Settings.SquareSize/Param.Settings.SquareLength;      

%*************************************
%*********Keys for Response***********
%*************************************
Param.Keys.Space        = 32;  
Param.Keys.EscPress     = 27;
Param.Keys.Left         = 37; 
Param.Keys.Right        = 39;
Param.Keys.one1         = 49; 
Param.Keys.two1         = 50;
Param.Keys.one2         = 97; 
Param.Keys.two2         = 98;
Param.Keys.Trigger1     = 83;  % 's'    


%*************************************
%*********Stimulus Locations**********
%*************************************
Param.Stimuli.Eccentricity   = 9; % degree

% location(1) = Left Visual Field;  location(2) = Right Visual Field;   location(3) = Center Visual Field;
Param.Stimuli.Locations(1,:) = [Param.Settings.ScrnResolution(3)/2 - Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree,Param.Settings.ScrnResolution(4)/2];
Param.Stimuli.Locations(2,:) = [Param.Settings.ScrnResolution(3)/2 + Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree,Param.Settings.ScrnResolution(4)/2];
Param.Stimuli.Locations(3,:) = [Param.Settings.ScrnResolution(3)/2, Param.Settings.ScrnResolution(4)/2];
Param.Stimuli.LocationsText  = {'Left','Right','Center'};


%*************************************
%******Parameters for Gratings********
%*************************************
%Param.Stimuli.OuterRadius      = 4;
%Param.Stimuli.InnerRadius      = Param.Stimuli.OuterRadius - 0.5;
%Param.Stimuli.GratingSize      = round(Param.Stimuli.OuterRadius*Param.Settings.PixelPerDegree);   % radius pixels
%Param.Stimuli.GratingOri       = 55;  
%Param.Stimuli.GratingNum       = size(Param.Stimuli.Locations,1); 
%Param.Stimuli.GratingContrast  = 0.8;
%Param.Stimuli.Spatial_freq     = 1/Param.Settings.PixelPerDegree;          % 0.02
%Param.Stimuli.SmoothSD         = Param.Settings.PixelPerDegree/6;          % 0.5 degree*1/3


%*************************************
%******Parameters for Fixation********
%*************************************
Param.Fixation.CrossColor    = [1,1,1]*255;
Param.Fixation.CrossSize     = 0.4*Param.Settings.PixelPerDegree;
Param.Fixation.CrossWidth    = 0.1*Param.Settings.PixelPerDegree;
Param.Fixation.CrossLoc      = [Param.Settings.ScrnResolution(3)/2-Param.Fixation.CrossSize, Param.Settings.ScrnResolution(3)/2+Param.Fixation.CrossSize, Param.Settings.ScrnResolution(3)/2, Param.Settings.ScrnResolution(3)/2;...
                                     Param.Settings.ScrnResolution(4)/2, Param.Settings.ScrnResolution(4)/2, Param.Settings.ScrnResolution(4)/2-Param.Fixation.CrossSize, Param.Settings.ScrnResolution(4)/2+Param.Fixation.CrossSize]; 

Param.Fixation.OvalSize      = 0.5*Param.Settings.PixelPerDegree;
Param.Fixation.OvalLoc       = [Param.Settings.ScrnResolution(3)/2-Param.Fixation.OvalSize, Param.Settings.ScrnResolution(4)/2-Param.Fixation.OvalSize,Param.Settings.ScrnResolution(3)/2+Param.Fixation.OvalSize,Param.Settings.ScrnResolution(4)/2+Param.Fixation.OvalSize];

Param.Fixation.OvalColor     = [0,0,0]; 
% Param.Fixation.ChangeProb    = 0.8;


%*************************************
%*********Parameters for RDK**********
%*************************************
Param.RDK.DotNum     = 400;
Param.RDK.DotSize    = 3; 

%Param.RDK.Direction  = 360-120;   % 0=Rihgt 90=down 180=left 270=up
Param.RDK.Speed      = 14.2;         % Watanabe 14.2  %点运动速度，用于控制实验难度

Param.RDK.OuterRadius= 5*Param.Settings.PixelPerDegree;         % radius
Param.RDK.InnerRadius= 0.5*Param.Settings.PixelPerDegree;       % radius
Param.RDK.Duration   = 0.4;                                     % s
Param.RDK.TesDuration   = 1;                                    % s
Param.RDK.FramesPerMove = 1;                                    % 每XX个frame更新点的位置
Param.RDK.NumFrames  = Param.RDK.Duration*NominalFrameRate/Param.RDK.FramesPerMove;         %注意：这里必须是整数
Param.RDK.StepPerMove   = Param.RDK.FramesPerMove/NominalFrameRate*Param.RDK.Speed*Param.Settings.PixelPerDegree; %移动pixel数
Param.RDK.TestNumFrames  = Param.RDK.TesDuration*NominalFrameRate/Param.RDK.FramesPerMove;  %注意：这里必须是整数

%*************************************
%*********Parameters for cue**********
%*************************************
Param.cue.Duration = 0.4;
Param.cue.Display = ['1','2'];

%*************************************
%*********Parameters for delay**********
%*************************************
Param.delay.Duration = 10;

%*************************************
%*******Discrimination Practice*******
%*************************************
Param.DisPrac.TrialNum     = 20;
Param.DisPrac.AngleDelta   = 5;          %discrimination
Param.DisPrac.Coherence    = 0.4;         %0-1

Param.DisPrac.Directions   = [15,75,135,195,255,315];% 6 directions
Param.DisPrac.DirectionNum = 6;
% Param.DisPrac.TrialPercond = 32;

Param.DisPrac.StiLocation  = 3;
Param.VarDir = 15;


%*************************************
%*******Discrimination Testing********
%*************************************
Param.DisTest.AngleDelta   = 5;          %discrimination
Param.DisTest.Coherence    = 1;          %0-1
Param.DisTest.StiLocation  = 3;


%*************************************
%**********Detection Practice*********
%*************************************
Param.DetectPrac.TrialNum     = 30;
Param.DetectPrac.Coherence    = 0.25;           %0-1

Param.DetectPrac.Directions   = [90,225,315];%3 directions
Param.DetectPrac.DirectionNum = 3;
Param.DetectPrac.TrialPercond = 10;

Param.DetectPrac.StiLocation  = 3;


%*************************************
%**********Detection Testing**********
%*************************************
Param.DetectTest.StiLocation  = 3;


%*************************************
%*********Detection Training**********
%*************************************
Param.DetectTrain.StiLocation  = 3;


%*************************************
%*******Parameters for Staircase******
%*************************************
Param.Staircase.AngleDelta          = 5;        %discrimination
Param.Staircase.Coherence           = 0.25;     %detection
Param.Staircase.MaxTrial            = 100;

Param.Staircase.Up                  = 1;        %increase after 1 wrong
Param.Staircase.Down                = 3;        %decrease after 3 consecutive right
Param.Staircase.StepSizeDown        = 0.25;     %discrimination        
Param.Staircase.StepSizeUp          = 0.25;     %discrimination
Param.Staircase.StepSizeDownDetect  = 0.01;    %detection        
Param.Staircase.StepSizeUpDetect    = 0.01;    %detection 

Param.Staircase.StopCriterion = 'reversals';   
Param.Staircase.StopRule      = 15;
Param.Staircase.ReversalsUsed = 7; 

Param.Staircase.xMax        = 15;               %discrimination
Param.Staircase.xMaxDetect  = 1;                %detection
Param.Staircase.xMin        = 0;


%*************************************
%*********Parameters for Trials*******
%*************************************
Param.Trial.ITI = 2; 
Param.Trial.ISI = 0.4; 
Param.Trial.MaxRT = 3;
% Param.Trial.Duration = 3;  % 0.6+0.2+0.5+0.2+1.4 = 2.9s
Param.Trial.Duration = 16.7;  % 0.4+0.4+0.4+0.4+0.4+10+1+3+2= 18s


%%%%%%%%%%%%%%%%%%%%提示声音%%%%%%%%%%%%%%%%%%%
ff=1500;%%声音频率
duration=100 ;%%声音持续的时间，单位毫秒
samplefrq=10000;%%声音信号的采样频率
nn=duration/1000;
%errorsound_data=cos(linspace(0,ff*2*pi*nn,samplefrq*nn)); %错误提示
errorsound_data=MakeBeep(ff, duration/1000, samplefrq);
errorsound_player=audioplayer(errorsound_data,samplefrq);%%声音的handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

