% sample
% Snd ('Play',beep3);
% signal pattern
[x,y]=meshgrid(-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize,-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize);
phase = rand*2*pi;
angle = rand*pi;
noise_angle = angle;
Noise_sti = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle)+x*cos(angle))+phase));
Noise_sti = Matrix_shuffle(Noise_sti);
Noise_sti(sqrt(x.^2 + y.^2) > Param.Stimuli.OuterSize) = gray;

[x,y]=meshgrid(-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize,-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize);
phase = rand*2*pi;
angle = 135/180*pi;
Signal_sti = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle)+x*cos(angle))+phase).*exp(-(x.^2 + y.^2)/2/(Param.Stimuli.SmoothSD.^2)));
Signal_sti(sqrt(x.^2 + y.^2) > Param.Stimuli.OuterSize) = gray;

% selection
overlay_temp = ones(size(x));
overlay_temp(sqrt(x.^2 + y.^2) > Param.Stimuli.OuterSize) = 0;
overlay_temp(sqrt(x.^2 + y.^2) < Param.Stimuli.InnerSize) = 0;
selection_temp = find(overlay_temp==1);

overlay_temp = zeros(size(x));

signal_pixels_num = round(sample_ratio* length(selection_temp));

rand_temp = randperm(length(selection_temp));
signal_pixels = rand_temp(1:signal_pixels_num);
overlay_temp(selection_temp(signal_pixels)) = 1;

Sti_final = Noise_sti.*(1-overlay_temp)+Signal_sti.*overlay_temp;


mm = Screen('MakeTexture', wnd, Sti_final);
Screen('DrawTexture', wnd, mm,[],[Param.Stimuli.Locations(location_used,1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(location_used,2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(location_used,1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(location_used,2)+Param.Stimuli.OuterSize]);
Screen('FillOval',wnd, black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
vbl = Screen('Flip',wnd);  %-Param.Screen.Slack


%% close all
is_true = 0;
while ~is_true
    [ifkey,RT_time,keyCode] = KbCheck;
    if keyCode(Param.Keys.Esc)
        is_true = 1;
    end
end

Screen('CloseAll');
reset_test_gamma;