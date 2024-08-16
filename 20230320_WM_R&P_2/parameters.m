% BY Ke Jia: Version 1_20220623 15:31
% list all the parameters used in this experiment

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
Param.Stimuli.Eccentricity   = 5; % degree
Param.Stimuli.Locations(1,:) = [Param.Settings.ScrnResolution(3)/2 - Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree,Param.Settings.ScrnResolution(4)/2];
Param.Stimuli.Locations(2,:) = [Param.Settings.ScrnResolution(3)/2 + Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree,Param.Settings.ScrnResolution(4)/2];
Param.Stimuli.Locations(3,:) = [Param.Settings.ScrnResolution(3)/2 + 0.8*Param.Stimuli.Eccentricity*Param.Settings.PixelPerDegree, Param.Settings.ScrnResolution(4)/2 + 0.6*Param.Stimuli.Eccentricity*Param.Settings.PixelPerDegree];
Param.Stimuli.Locations(4,:) = [Param.Settings.ScrnResolution(3)/2 - 0.8*Param.Stimuli.Eccentricity*Param.Settings.PixelPerDegree, Param.Settings.ScrnResolution(4)/2 + 0.6*Param.Stimuli.Eccentricity*Param.Settings.PixelPerDegree];
Param.Stimuli.Locations(5,:) = [Param.Settings.ScrnResolution(3)/2 , Param.Settings.ScrnResolution(4)/2- Param.Stimuli.Eccentricity*Param.Settings.PixelPerDegree];
Param.Stimuli.Locations(6,:) = [Param.Settings.ScrnResolution(3)/2 , Param.Settings.ScrnResolution(4)/2];
Param.Stimuli.LocationsText  = {'Left','Right','Lower Right','Lower Left','Up','Centre'};
Param.Stimuli.Location_used  = [3 4 5];

% Gratings
Param.Stimuli.OuterRadius      = 2.5;
% Param.Stimuli.OuterRadius_edge = Param.Stimuli.OuterRadius - 0.5;
Param.Stimuli.OuterSize        = round(Param.Stimuli.OuterRadius*Param.Settings.PixelPerDegree);   % radius pixels
% Param.Stimuli.OuterSize2       = round(Param.Stimuli.OuterRadius2*Param.Settings.PixelPerDegree);   % radius pixels

Param.Stimuli.GratingOri       = [15 75 135];
Param.Stimuli.Orireward        = Param.Stimuli.GratingOri(Ori_reward);
Param.Stimuli.Oripunish        = Param.Stimuli.GratingOri(length(Param.Stimuli.GratingOri) - Ori_reward);
Param.Stimuli.OriJitter        = 10; 
Param.Stimuli.ResponseJitter   = 10; 

Param.Stimuli.GratingContrast  = 0.8;
Param.Stimuli.Spatial_freq     = 0.6/Param.Settings.PixelPerDegree;       % 0.02 % 6 bars in total
% Param.Stimuli.SmoothSD         = Param.Stimuli.OuterRadius*Param.Settings.PixelPerDegree/3;          % 0.5degree*1/3
 
%% parameters for fixation
Param.Fixation.CrossColor    = [1,1,1]*255;
Param.Fixation.CrossSize     = 0.3*Param.Settings.PixelPerDegree;
Param.Fixation.CrossWidth    = 0.1*Param.Settings.PixelPerDegree;
% +
Param.Fixation.CrossLoc      = [Param.Settings.ScrnResolution(3)/2-Param.Fixation.CrossSize, Param.Settings.ScrnResolution(3)/2+Param.Fixation.CrossSize, Param.Settings.ScrnResolution(3)/2, Param.Settings.ScrnResolution(3)/2;...
                                     Param.Settings.ScrnResolution(4)/2, Param.Settings.ScrnResolution(4)/2, Param.Settings.ScrnResolution(4)/2-Param.Fixation.CrossSize, Param.Settings.ScrnResolution(4)/2+Param.Fixation.CrossSize] + [offset',offset',offset',offset'];
% X
Param.Fixation.CrossLoc2     = [Param.Settings.ScrnResolution(3)/2-Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(3)/2+Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(3)/2-Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(3)/2+Param.Fixation.CrossSize*sqrt(2)/2;...
                                     Param.Settings.ScrnResolution(4)/2-Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(4)/2+Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(4)/2+Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(4)/2-Param.Fixation.CrossSize*sqrt(2)/2]+[offset',offset',offset',offset'];  

Param.Fixation.OvalSize      = 1*Param.Settings.PixelPerDegree;
Param.Fixation.OvalColor     = [0,0,0]; 
Param.Fixation.Offset        = 0.5*Param.Settings.PixelPerDegree; 
Param.Fixation.FontColor     = [1,1,1]*255;

%% parameters for trials
Param.Trial.ITI              = 0.4;     
Param.Trial.StiDura          = 2; 
Param.Trial.Delay            = 2;

Param.Trial.TestNum          = 72; 
Param.Trial.Testminirun      = 36;
Param.Trial.TrainNum         = 72; 
Param.Trial.Trainminirun     = 24;
Param.Trial.Practice         = 36;

Param.Trial.Response         = 1;  % 

%% open window 
screens = Screen('Screens');
screenNumber = max(screens);	
Screen('Preference', 'SkipSyncTests', 1);

white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
gray = round((white+black)/2);
wnd = Screen('OpenWindow',screenNumber,gray,Param.Settings.ScrnResolution); 
Screen('BlendFunction',wnd, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('TextFont',wnd,'Arial');
Screen('TextSize',wnd,30);

RefreshDur = Screen('GetFlipInterval',wnd);
RefreshRate = 1./RefreshDur;
Slack = RefreshDur/2;

% if abs(RefreshRate-100)>1
%     disp(' ');
%     disp('Please reset your refreshrate!');
%     disp(' ');
%     abort;
% end
