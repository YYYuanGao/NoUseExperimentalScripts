% BY YuanGao.
% list all the parameters used in this experiment.
%% 
Param = struct;

%% Screen Settings
Param.Settings.ViewDistance      = 750;              % 1100 mm 
Param.Settings.ScrnResolution    = [0 0 1024 768];   % rect_computer = [0 0 40 30];
Param.Settings.SquareSize        = Param.Settings.ScrnResolution(3);  % 400 mm
Param.Settings.SquareLength      = 360;               % 154 mm 
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
Param.Stimuli.Eccentricity   = 6; % degree
Param.Stimuli.Locations(1,:) = [Param.Settings.ScrnResolution(3)/2-Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree,Param.Settings.ScrnResolution(4)/2];
Param.Stimuli.Locations(2,:) = [Param.Settings.ScrnResolution(3)/2 + Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree,Param.Settings.ScrnResolution(4)/2];
Param.Stimuli.Locations(3,:) = [Param.Settings.ScrnResolution(3)/2 - 0.8*Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree, Param.Settings.ScrnResolution(4)/2+0.6*Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree];
Param.Stimuli.Locations(4,:) = [Param.Settings.ScrnResolution(3)/2 + 0.8*Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree, Param.Settings.ScrnResolution(4)/2+0.6*Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree];
Param.Stimuli.Locations(5,:) = [Param.Settings.ScrnResolution(3)/2 , Param.Settings.ScrnResolution(4)/2];
Param.Stimuli.LocationsText  = {'Left','Right','Lower Left','Lower Right','Centre'};
Param.Stimuli.Offset         = 10;
Param.Stimuli.VerOffset      = [-0.5 -0.25 0.25 0.5];
%% parameters for fixation
Param.Fixation.CrossColor    = [1,1,1]*255;
Param.Fixation.JitterColor   = [255,255,0];
Param.Fixation.CrossSize     = 0.3*Param.Settings.PixelPerDegree;
Param.Fixation.CrossWidth    = 0.1*Param.Settings.PixelPerDegree;
% +
Param.Fixation.CrossLoc      = [Param.Settings.ScrnResolution(3)/2-Param.Fixation.CrossSize, Param.Settings.ScrnResolution(3)/2+Param.Fixation.CrossSize, Param.Settings.ScrnResolution(3)/2, Param.Settings.ScrnResolution(3)/2;...
                                     Param.Settings.ScrnResolution(4)/2, Param.Settings.ScrnResolution(4)/2, Param.Settings.ScrnResolution(4)/2-Param.Fixation.CrossSize, Param.Settings.ScrnResolution(4)/2+Param.Fixation.CrossSize];
Param.Fixation.CrossLoc2     = [Param.Settings.ScrnResolution(3)/2-Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(3)/2+Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(3)/2-Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(3)/2+Param.Fixation.CrossSize*sqrt(2)/2;...
                                     Param.Settings.ScrnResolution(4)/2-Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(4)/2+Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(4)/2+Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(4)/2-Param.Fixation.CrossSize*sqrt(2)/2];  
Param.Fixation.OvalSize      = 1*Param.Settings.PixelPerDegree;
Param.Fixation.OvalColor     = [0,0,0]; 

%% Vernier
% Gratings
Param.Stimuli.OuterRadius      = 3;
% Param.Stimuli.OuterRadius_edge = Param.Stimuli.OuterRadius - 0.5;
Param.Stimuli.InnerRadius      = 0.3;
Param.Stimuli.OuterSize        = round(Param.Stimuli.OuterRadius*Param.Settings.PixelPerDegree);   % radius pixels
Param.Stimuli.InnerSize        = round(Param.Stimuli.InnerRadius*Param.Settings.PixelPerDegree);   % radius pixels

Param.Stimuli.GratingOri       = 0;
Param.Stimuli.OriJitter        = 5; 

Param.Stimuli.GratingContrast  = 0.7;
Param.Stimuli.Spatial_freq     = 5/Param.Settings.PixelPerDegree;       % 0.02 % 6 bars in total
Param.Stimuli.SmoothSD         = Param.Settings.PixelPerDegree/3;       % 0.5degree*1/3 Param.Stimuli.OuterRadius*

%% parameters for trials
Param.Trial.ITI               = 0.5;     
Param.Trial.StiDura1          = 2; 
Param.Trial.StiDura2          = [0.8 1.2 1.5]; 
Param.Trial.MaxRT             = 1.5; 
Param.Trial.Duration          = 6;
Param.Trial.RunsNum           = 300;
Param.Trial.MinNum            = 10;

Param.Trial.TestITI           = [0.3 0.4 0.5];     
Param.Trial.MaxRT             = 1;
Param.Trial.TestStiDura1      = 0.5; 
Param.Trial.TestStiDura2      = 0.3;
Param.Trial.TestMask          = 0.2;
Param.Trial.TestDuration      = 3;
Param.Trial.TestRunsNum       = 450;
Param.Trial.TestMinNum        = 150;
%% open window 
screens = Screen('Screens');
screenNumber = max(screens);	
Screen('Preference', 'SkipSyncTests', 0);

white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
gray = round((white+black)/2);
wnd = Screen('OpenWindow',screenNumber,gray,Param.Settings.ScrnResolution); 
Screen('BlendFunction',wnd,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

RefreshDur = Screen('GetFlipInterval',wnd);
Slack = RefreshDur/2;

