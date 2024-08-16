% list all the parameters used in this experiment
%% 
Param = struct;
% sp = BioSemiSerialPort();
%% Marker
Param.start                = 77;
Param.ITI1                 = 10;
Param.sti1                 = [101 102 103 104 105 106];
Param.ISI                  = 30;
Param.sti2                 = 200;
Param.response             = 40;
Param.correct              = 51;
Param.wrong                = 52;
Param.ITI2                 = 60;
%% Screen Settings
Param.Settings.ViewDistance      = 900;              % 1100 mm 
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
Param.Stimuli.Eccentricity   = 6; % degree
Param.Stimuli.Locations(1,:) = [Param.Settings.ScrnResolution(3)/2 - Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree,Param.Settings.ScrnResolution(4)/2];
Param.Stimuli.Locations(2,:) = [Param.Settings.ScrnResolution(3)/2 + Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree,Param.Settings.ScrnResolution(4)/2];
Param.Stimuli.Locations(3,:) = [Param.Settings.ScrnResolution(3)/2 + 0.8*Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree, Param.Settings.ScrnResolution(4)/2+0.6*Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree];
Param.Stimuli.Locations(4,:) = [Param.Settings.ScrnResolution(3)/2 , Param.Settings.ScrnResolution(4)/2];
Param.Stimuli.LocationsText  = {'Left','Right','Lower Right','Centre'};

% Gratings
Param.Stimuli.OuterRadius      = 5;
% Param.Stimuli.OuterRadius_edge = Param.Stimuli.OuterRadius - 0.5;
Param.Stimuli.InnerRadius      = 4.5;
Param.Stimuli.OuterSize        = round(Param.Stimuli.OuterRadius*Param.Settings.PixelPerDegree);   % radius pixels
Param.Stimuli.InnerSize        = round(Param.Stimuli.InnerRadius*Param.Settings.PixelPerDegree);   % radius pixels

Param.Stimuli.GratingOri       = [15 75 135 195 255 315];
Param.Stimuli.OriJitter        = 5; 

Param.Stimuli.GratingContrast  = 0.8;
Param.Stimuli.Spatial_freq     = 0.75/Param.Settings.PixelPerDegree;       % 0.02 % 6 bars in total
Param.Stimuli.SmoothSD         = Param.Settings.PixelPerDegree/6;       % 0.5degree*1/3 Param.Stimuli.OuterRadius*
 
%% parameters for testing_Quest
Param.RF.alphas = 0:0.5:15;  %-3:.01:3
Param.RF.sd = 2;
Param.RF.prior = PAL_pdfNormal(Param.RF.alphas,5,Param.RF.sd);
Param.RF.stopcriterion = 'trials';
Param.RF.stoprule = 10;  % 10 runs
Param.RF.stoprule1 = 5;  % 2 runs
Param.RF.PFfit = @PAL_Gumbel;    %Shape to be assumed
Param.RF.gamma = 0.5;               
Param.RF.beta = 3.5;               %Slope to be assumed
Param.RF.lambda  = 0.02;            %Lapse rate to be assumed
Param.RF.meanmode = 'mean';      %Use mean of posterior as placement rule
Param.RF.start = 5;
Param.RF.xMax = 15;
Param.RF.xMin = 0;
Param.RF.trueParams = [Param.RF.alphas Param.RF.beta Param.RF.gamma Param.RF.lambda];
Param.RF.PFsimul = @PAL_Gumbel;

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
Param.Trial.ITI1             = [0.5 0.75 1];
Param.Trial.ITI2             = 0.5;     
Param.Trial.ISI              = 1.3;    
Param.Trial.StiDura          = 0.2; 
Param.Trial.MaxRT            = 1.5; 
Param.Trial.Duration         = 5;
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
