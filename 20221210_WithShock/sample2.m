clc;
clear;
HideCursor;
parameters;
%% load
load pic;
for a = 1:length(pic)
    pic{a,3} = imread(pic{a,2});
    pic{a,3}(pic{a,3} < 49) = 128;
%     pic{a,3}(pic{a,3} == 34) = 128;
%     pic{a,3}(pic{a,3} == 0) = 128;
%     pic{a,3}(pic{a,3} == 0) = 128;
%     pic{a,3}(pic{a,3} == 0) = 128;
%     pic{a,3}(pic{a,3} == 0) = 128;
%     pic{a,3}(pic{a,3} == 0) = 128;
%     pic{a,3}(pic{a,3} == 0) = 128;
%     pic{a,3}(pic{a,3} == 0) = 128;
%     pic{a,3}(pic{a,3} == 0) = 128;
%     pic{a,3}(pic{a,3} == 0) = 128;
%     pic{a,3}(pic{a,3} == 0) = 128;
    
end

%% maketex
i = randperm(length(pic),1);
j = randperm(length(pic),1);
stim1 = Screen('MakeTexture',wnd,pic{i,3});  
stim2 = Screen('MakeTexture',wnd,pic{j,3});
%% loc
width = size(pic{a,3},2);
height = size(pic{a,3},1);
stim1_loc = [Param.Settings.ScrnResolution(3)/4-width/2,Param.Settings.ScrnResolution(4)/2-height/2,Param.Settings.ScrnResolution(3)/4+width/2,Param.Settings.ScrnResolution(4)/2+height/2];
stim2_loc = [3*Param.Settings.ScrnResolution(3)/4-width/2,Param.Settings.ScrnResolution(4)/2-height/2,3*Param.Settings.ScrnResolution(3)/4+width/2,Param.Settings.ScrnResolution(4)/2+height/2];

%% draw
location_temp1 = [1 2];
loc1 = randperm(length(location_temp1),1);
location_temp1 = location_temp1(loc1);

if location_temp1 == 1
    location_temp2 = 2;
elseif location_temp1 == 2
    location_temp2 = 1;
    loc2 = randperm(length(location_temp2),1);
    location_temp2 = location_temp2(loc2);
end
Screen('DrawTexture', wnd, stim1, [],[Param.Stimuli.Locations(location_temp1,1)-width/3,Param.Stimuli.Locations(location_temp1,2)-height/3,Param.Stimuli.Locations(location_temp1,1)+width/3,Param.Stimuli.Locations(location_temp1,2)+height/3]);
Screen('DrawTexture', wnd, stim2, [],[Param.Stimuli.Locations(location_temp2,1)-width/3,Param.Stimuli.Locations(location_temp2,2)-height/3,Param.Stimuli.Locations(location_temp2,1)+width/3,Param.Stimuli.Locations(location_temp2,2)+height/3]);

% fixation
Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, white, [], 1);
Screen('Flip',wnd);
WaitSecs(0.5);
%% noise
t = imread('gray.png');
t1 = imnoise(t,'gaussian',0,0.08);
save t1
load t1
noise_stim = Screen('MakeTexture',wnd,t1);  
Screen('DrawTexture',wnd,noise_stim);                            
Screen('Flip',wnd);
WaitSecs(0.15);

%% 
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
%% 
[x,y]=meshgrid(-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize,-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize);
phase = rand*2*pi;
angle = Param.Stimuli.GratingOri/180*pi;
annulus=ones(Param.Stimuli.OuterSize*2+1);
edge_control = sqrt(x.^2 + y.^2)./Param.Settings.PixelPerDegree;
overrado=find(edge_control>Param.Stimuli.InnerRadius);

len=length(overrado);
for i=1:len
    annulus(overrado(i))=(annulus(overrado(i)).*exp(-1.*((((edge_control(overrado(i))-Param.Stimuli.InnerRadius)*Param.Settings.PixelPerDegree).^2)/(2.*(Param.Stimuli.SmoothSD.^2)))));
end

% sin grating
Gabor_sti = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle)+x*cos(angle))+phase).*exp(-(x.^2 + y.^2)/(2.*Param.Stimuli.SmoothSD.^2)));
Gabor_sti(sqrt(x.^2 + y.^2) > Param.Stimuli.OuterSize) = gray;

Gabor_sti_final = repmat(Gabor_sti,[1,1,3]);
Gabor_sti_final(:,:,4) = annulus*white;

mm = Screen('MakeTexture', wnd, Gabor_sti_final);

location_temp3 = [5 10];
loc3 = randperm(length(location_temp3),1);
location_temp3 = location_temp3(loc3);

if location_temp3 == 5
    location_temp4 = [6 7 8 9];
    loc4 = randperm(length(location_temp4),1);
    location_temp4 = location_temp4(loc4);

elseif location_temp3 == 10
    location_temp4 = [11 12 13 14];
    loc4 = randperm(length(location_temp4),1);
    location_temp4 = location_temp4(loc4);
end



Screen('DrawTexture', wnd, mm,[],[Param.Stimuli.Locations(location_temp3,1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(location_temp3,2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(location_temp3,1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(location_temp3,2)+Param.Stimuli.OuterSize]);
Screen('DrawTexture', wnd, mm,[],[Param.Stimuli.Locations(location_temp4,1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(location_temp4,2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(location_temp4,1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(location_temp4,2)+Param.Stimuli.OuterSize]);

% fixation
Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, white, [], 1);
Screen('Flip',wnd); %-Slack

     
%% close all
is_true = 0;
while ~is_true
    [~,~,keyCode] = KbCheck;
    if keyCode(Param.Keys.Esc)
        is_true = 1;
    end
end

Screen('CloseAll');
reset_test_gamma;