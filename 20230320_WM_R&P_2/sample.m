[x,y]=meshgrid(-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize,-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize);
angle1 = Param.Stimuli.Orireward/180*pi;   
angle2 = Param.Stimuli.Oripunish/180*pi;   
angle3 = Param.Stimuli.GratingOri(3)/180*pi;  

phase = rand*2*pi;
sti1_final = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle1)+x*cos(angle1))+phase));
sti1_final(sqrt(x.^2 + y.^2) > Param.Stimuli.OuterSize) = gray; %circle mask

phase = rand*2*pi;
sti2_final = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle2)+x*cos(angle2))+phase));
sti2_final(sqrt(x.^2 + y.^2) > Param.Stimuli.OuterSize) = gray; %circle mask

phase = rand*2*pi;
sti3_final = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle3)+x*cos(angle3))+phase));
sti3_final(sqrt(x.^2 + y.^2) > Param.Stimuli.OuterSize) = gray; %circle mask

mm1 = Screen('MakeTexture', wnd, sti1_final);
mm2 = Screen('MakeTexture', wnd, sti2_final);
mm3 = Screen('MakeTexture', wnd, sti3_final);


location_temp1 = Param.Stimuli.Location_used(1);
location_temp2 = Param.Stimuli.Location_used(2);
location_temp3 = Param.Stimuli.Location_used(3);

Screen('DrawTexture', wnd, mm1,[],[Param.Stimuli.Locations(location_temp1,1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(location_temp1,2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(location_temp1,1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(location_temp1,2)+Param.Stimuli.OuterSize]);
Screen('DrawTexture', wnd, mm2,[],[Param.Stimuli.Locations(location_temp2,1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(location_temp2,2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(location_temp2,1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(location_temp2,2)+Param.Stimuli.OuterSize]);
Screen('DrawTexture', wnd, mm3,[],[Param.Stimuli.Locations(location_temp3,1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(location_temp3,2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(location_temp3,1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(location_temp3,2)+Param.Stimuli.OuterSize]);

% fixation
Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
vbl = Screen('Flip',wnd);

%% close all
is_true = 0;
while ~is_true
    [~,~,keyCode] = KbCheck;
    if keyCode(Param.Keys.Esc)
        is_true = 1;
    end
end

Screen('CloseAll');
