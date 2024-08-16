% list all the parameters used in this experiment
%% 
Param = struct;

%% Screen Settings
Param.Settings.ViewDistance    = 600;              % 1100 mm 
Param.Settings.ScrnResolution  = [0 0 1024 768];   % rect_computer = [0 0 40 30];
Param.Settings.SquareSize      = Param.Settings.ScrnResolution(3);             % 400 mm 
Param.Settings.SquareLength    = 360;              % 154 mm  
Param.Settings.PixelPerDegree  = 2*Param.Settings.ViewDistance*tan(1/2/180*pi)*Param.Settings.SquareSize/Param.Settings.SquareLength;       

%% Keys for response
Param.Keys.Space    = 32;  
Param.Keys.Esc      = 27;
Param.Keys.Left     = 37;  
Param.Keys.Right    = 39;
Param.Keys.Down     = 40;
Param.Keys.Trigger1 = 83;  % 's'       

%% parameters for stimuli
% Locations
Param.Stimuli.Eccentricity   = 6; % degree
Param.Stimuli.Locations(1,:) = [Param.Settings.ScrnResolution(3)/2 - Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree,Param.Settings.ScrnResolution(4)/2];
Param.Stimuli.Locations(2,:) = [Param.Settings.ScrnResolution(3)/2 + Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree,Param.Settings.ScrnResolution(4)/2];
Param.Stimuli.Locations(3,:) = [Param.Settings.ScrnResolution(3)/2 + 0.8*Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree, Param.Settings.ScrnResolution(4)/2+0.6*Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree];
Param.Stimuli.Locations(4,:) = [Param.Settings.ScrnResolution(3)/2 , Param.Settings.ScrnResolution(4)/2];
Param.Stimuli.LocationsText  = {'Left','Right','Lower Right','Centre'};

% Gratings
Param.Stimuli.OuterRadius      = 5;
Param.Stimuli.OuterRadius_edge = Param.Stimuli.OuterRadius - 0.5;
Param.Stimuli.InnerRadius      = 1;
Param.Stimuli.OuterSize        = round(Param.Stimuli.OuterRadius*Param.Settings.PixelPerDegree);   % radius pixels
Param.Stimuli.InnerSize        = round(Param.Stimuli.InnerRadius*Param.Settings.PixelPerDegree);   % radius pixels

Param.Stimuli.GratingOri       = [15 75 135];
Param.Stimuli.OriJitter        = 5; 

Param.Stimuli.GratingContrast  = 0.8;
Param.Stimuli.Spatial_freq     = 0.6/Param.Settings.PixelPerDegree;       % 0.02 % 6 bars in total
Param.Stimuli.SmoothSD         = Param.Stimuli.OuterRadius*Param.Settings.PixelPerDegree/2.3;       % 0.5degree*1/3 Param.Stimuli.OuterRadius*

%% parameters for cues
% InitializePsychSound;
% freq = 48000;
% pahandle = PsychPortAudio('Open', [], 1, 1, freq,2);
% PsychPortAudio('Volume', pahandle, 0.5);
% [beep0,samplingRate] = MakeBeep(1000,0.1,freq);
% [beep1,samplingRate] = MakeBeep(800,0.1,freq);
% [beep2,samplingRate] = MakeBeep(1000,0.1,freq);
% [beep3,samplingRate] = MakeBeep(1200,0.1,freq);

%% parameters for testing_staircase
Param.Stairtest.start = 5;
Param.Stairtest.MaxTrial = 40;
Param.Stairtest.Up = 1;     %increase after 1 wrong
Param.Stairtest.Down = 2;   %decrease after 2 consecutive right
Param.Stairtest.StepSizeDown = 0.25;         
Param.Stairtest.StepSizeUp = 0.25;
Param.Stairtest.StopCriterion = 'trials';   
Param.Stairtest.StopRule1 = 40;
Param.Stairtest.StopRule2 = 20;
Param.Stairtest.ReversalsUnused = 3; 
Param.Stairtest.xMax = 15;
Param.Stairtest.xMin = 0;
Param.Stairtest.trialthresh = 5;

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

%% parameters for trials
Param.Trial.ITI              = [0.5 1 1.5];     
Param.Trial.ISI              = 0.5;    
Param.Trial.StiDura          = 0.15; 
Param.Trial.MaxRT            = 1.3; 
Param.Trial.Duration         = 4; 
%% open window 
screens = Screen('Screens');
screenNumber = max(screens);	
Screen('Preference', 'SkipSyncTests', 0);

white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
gray = round((white+black)/2);
wnd = Screen('OpenWindow',screenNumber,gray,Param.Settings.ScrnResolution); 
Screen('BlendFunction',wnd, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% RefreshDur = Screen('GetFlipInterval',wnd);
% Slack = RefreshDur/2;
% RefreshRate = 1/RefreshDur;
% if abs(RefreshRate-100)>1
%     disp(' ');
%     disp('Please double check the RefreshRate!');
%     disp(' ');
%     abort;
% end
