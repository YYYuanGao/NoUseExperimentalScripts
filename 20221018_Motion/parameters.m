% list all the parameters used in this experiment
%% 
Param = struct;

%% Screen Settings
Param.Settings.ViewDistance      = 500;              % 1100 mm 
Param.Settings.ScrnResolution    = [0 0 1920 1080];   % rect_computer = [0 0 40 30];
Param.Settings.SquareSize        = Param.Settings.ScrnResolution(3);  % 400 mm
Param.Settings.SquareLength      = Screen('DisplaySize',0);               % 154 mm 
Param.Settings.PixelPerDegree    = 2*Param.Settings.ViewDistance*tan(1/2/180*pi)*Param.Settings.SquareSize/Param.Settings.SquareLength;       

%% Keys for response
Param.Keys.Space    = 32;  
Param.Keys.Esc      = 27;
Param.Keys.Left     = 37;  
Param.Keys.Right    = 39;
Param.Keys.Down     = 40;
Param.Keys.Trigger1 = 83;  % 's'       

%% parameters for stimuli
% Locations
Param.Stimuli.Eccentricity   = 9; % degree
Param.Stimuli.Locations(1,:) = [Param.Settings.ScrnResolution(3)/2 - Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree,Param.Settings.ScrnResolution(4)/2];
Param.Stimuli.Locations(2,:) = [Param.Settings.ScrnResolution(3)/2 + Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree,Param.Settings.ScrnResolution(4)/2];
Param.Stimuli.Locations(3,:) = [Param.Settings.ScrnResolution(3)/2 + 0.8*Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree, Param.Settings.ScrnResolution(4)/2+0.6*Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree];
Param.Stimuli.Locations(4,:) = [Param.Settings.ScrnResolution(3)/2 , Param.Settings.ScrnResolution(4)/2];
Param.Stimuli.LocationsText  = {'Left','Right','Lower Right','Centre'};

%% parameters for fixation
Param.Fixation.CrossColor    = [1,1,1]*255;
Param.Fixation.CrossSize     = 0.3*Param.Settings.PixelPerDegree;
Param.Fixation.CrossWidth    = 0.1*Param.Settings.PixelPerDegree;
% +
Param.Fixation.CrossLoc      = [Param.Settings.ScrnResolution(3)/2-Param.Fixation.CrossSize, Param.Settings.ScrnResolution(3)/2+Param.Fixation.CrossSize, Param.Settings.ScrnResolution(3)/2, Param.Settings.ScrnResolution(3)/2;...
                                     Param.Settings.ScrnResolution(4)/2, Param.Settings.ScrnResolution(4)/2, Param.Settings.ScrnResolution(4)/2-Param.Fixation.CrossSize, Param.Settings.ScrnResolution(4)/2+Param.Fixation.CrossSize];
% X
Param.Fixation.CrossLoc2     = [Param.Settings.ScrnResolution(3)/2-Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(3)/2+Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(3)/2-Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(3)/2+Param.Fixation.CrossSize*sqrt(2)/2;...
                                     Param.Settings.ScrnResolution(4)/2-Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(4)/2+Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(4)/2+Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(4)/2-Param.Fixation.CrossSize*sqrt(2)/2];  
Param.Fixation.OvalSize      = 1*Param.Settings.PixelPerDegree;
Param.Fixation.OvalColor     = [0,0,0]; 

%% Staircase
Param.Staircase.AngleDelta          = 5;        
Param.Staircase.MaxTrial            = 100;
Param.Staircase.Up                  = 1;        %increase after 1 wrong
Param.Staircase.Down                = 3;        %decrease after 3 consecutive right
Param.Staircase.StepSizeDown        = 0.25;             
Param.Staircase.StepSizeUp          = 0.25;     
Param.Staircase.StopCriterion       = 'reversals';   
Param.Staircase.StopRule            = 15;
Param.Staircase.ReversalsUsed       = 7; 
Param.Staircase.xMax                = 15;             
Param.Staircase.xMin                = 0;

%% Parameters for RDK
Param.RDK.DotNum                    = 100;
Param.RDK.DotSize                   = 3; 

% Param.RDK.Direction                 = 360-120;                                 % 0=Rihgt 90=down 180=left 270=up
Param.RDK.Speed                     = 14.2;         
Param.RDK.OuterRadius               = 5*Param.Settings.PixelPerDegree;         % radius
Param.RDK.InnerRadius               = 0.5*Param.Settings.PixelPerDegree;       % radius
Param.RDK.Duration                  = 0.4;                                     % s
Param.RDK.TesDuration               = 1;                                       % s
Param.RDK.FramesPerMove             = 1;                                       % 每XX个frame更新点的位置
Param.RDK.NumFrames                 = Param.RDK.Duration*NominalFrameRate/Param.RDK.FramesPerMove;                                %Integer
Param.RDK.StepPerMove               = Param.RDK.FramesPerMove/NominalFrameRate*Param.RDK.Speed*Param.Settings.PixelPerDegree;     %moving pixel
Param.RDK.TestNumFrames             = Param.RDK.TesDuration*NominalFrameRate/Param.RDK.FramesPerMove;                             %Integer


%% Test
Param.DisTest.Coherence                   = 1;          %0-1
% Param.DisPrac.TrialNum                  = 20;
Param.Trial.Directions                    = [15,75,135];
% Param.DisPrac.DirectionNum              = 6;
% % Param.DisPrac.TrialPercond            = 32;
% Param.VarDir                            = 15;
%% parameters for trials
Param.Trial.ITI                         = 2;     
Param.Trial.ISI                         = 0.5;    
Param.Trial.StiDura                     = 0.15; 
Param.Trial.MaxRT                       = 1.3; 
Param.Trial.Duration                    = 16.7;     % 0.4+0.4+0.4+0.4+0.4+10+1+3+2= 18s
%% Open Window
screens = Screen('Screens');
screenNumber = max(screens);	

white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
gray = round((white+black)/2);
Screen('Preference', 'SkipSyncTests', 0);
wnd = Screen('OpenWindow',screenNumber,gray,Param.Settings.ScrnResolution); 
Screen('BlendFunction', wnd, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

FrameRate = Screen('FrameRate',wnd);
NominalFrameRate = Screen('NominalFrameRate',wnd); 
RefreshDur = Screen('GetFlipInterval',wnd);
Slack = RefreshDur / 2;

% RefreshDur = Screen('GetFlipInterval',wnd);
% Slack = RefreshDur/2;
% RefreshRate = 1/RefreshDur;
% if abs(RefreshRate-100)>1
%     disp(' ');
%     disp('Please double check the RefreshRate!');
%     disp(' ');
%     abort;
% end

