parameters;

    %% display sample1
    for frame_i = 1:Param.RDK.NumFrames
        %         Screen('DrawDots',wnd,MyDot_1{1,frame_i},Param.RDK.DotSize,black,Param.Stimuli.Locations(1,:),1);
        Screen('DrawDots',wnd,MyDot_1{1,frame_i},Param.RDK.DotSize,black,Param.Stimuli.Locations(stimLoc,:),1);
        Screen('FillOval', wnd, Param.Fixation.OvalColor, Param.Fixation.OvalLoc);
        Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);


        if frame_i == 1
            vbl = Screen('Flip',wnd,vbl+Param.Trial.ITI-Slack);
            results(trial_i,14)=vbl;  %%motion onset
        else
            vbl = Screen('Flip',wnd,vbl+Param.RDK.FramesPerMove/NominalFrameRate-Slack);
        end
    end
    
    %% ISI
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
    vbl = Screen('Flip',wnd,vbl+Param.Trial.StiDura); %-Slack
    results(trial_i,14) = vbl-trial_onset-sum(results(trial_i,12:13)); 

    %% sti2
    [x,y]=meshgrid(-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize,-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize);
    phase = rand*2*pi; 
    angle = results(trial_i,5)/180*pi;
    annulus=ones(Param.Stimuli.OuterSize*2+1);
    edge_control = sqrt(x.^2 + y.^2)./Param.Settings.PixelPerDegree;
    overrado=find(edge_control>Param.Stimuli.InnerRadius);

    len=length(overrado);
    for i=1:len
        annulus(overrado(i))=(annulus(overrado(i)).*exp(-1.*((((edge_control(overrado(i))-Param.Stimuli.InnerRadius)*Param.Settings.PixelPerDegree).^2)/((2.*Param.Stimuli.SmoothSD.^2)))));
    end













location_temp = 4;
Screen('DrawTexture', wnd, mm,[],[Param.Stimuli.Locations(location_temp,1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(location_temp,2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(location_temp,1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(location_temp,2)+Param.Stimuli.OuterSize]);


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
reset_test_gamma;