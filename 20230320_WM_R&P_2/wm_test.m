% WM with Reward and Punishment.
% By YuanGao. 20230314 18:12

%% check parameters used
if exist([CurrDir '\Results\test\' SubjID '\' SubjID '_results_sess' num2str(sess_num) '_run' num2str(run_num) '.mat'],'file')
    disp(' ');
    disp([SubjID '_run' num2str(run_num) ' has been test, please enter a new run number!']);
    disp(' ');
    abort;
end

%% 
results = zeros(Param.Trial.TestNum,22);
trial_index = randperm(Param.Trial.TestNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% results:
% 1 trial number            % 10 ori for control_show       % 19 ori error 
% 2 loc for reward          % 11 cued_condition             % 20 trial duration 
% 3 loc for punish          % 12 reported ori               % 21 RT
% 4 loc for control         % 13 trial onset                % 22 cued_loc      
% 5 ori for reward          % 14 fix onset delay            
% 6 ori for punish          % 15 ITI
% 7 ori for control         % 16 stim
% 8 ori for reward_show     % 17 delay
% 9 ori for punish_show     % 18 report startpoint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Go!
for trial_i = 1:Param.Trial.TestNum
    
    % start
    if mod(trial_i,Param.Trial.Testminirun) == 1
        if trial_i == 1
            DrawFormattedText(wnd,'Press space to start!','center','center', black);
        else
            DrawFormattedText(wnd,'Take a rest! Press space to continue!' ,'center','center',black);
            WaitSecs(2);
        end
    
        Screen('Flip',wnd);
        is_true = 0;
        while (is_true == 0)
            [ifkey,RT_time,keyCode] = KbCheck;
            if keyCode(Param.Keys.Space)
                is_true = 1;
            elseif keyCode(Param.Keys.Esc)
                abort;          
            end
        end
    end

    results(trial_i,1) = trial_i;
    temp1 = randperm(length(Param.Stimuli.Location_used));
    results(trial_i,2) = Param.Stimuli.Location_used(temp1(1));  % loc for reward
    results(trial_i,3) = Param.Stimuli.Location_used(temp1(2));  % loc for punish
    results(trial_i,4) = Param.Stimuli.Location_used(temp1(3));  % loc for control

    % define ref orientation
    results(trial_i,5) = Param.Stimuli.Orireward;            % ori for reward
    results(trial_i,6) = Param.Stimuli.Oripunish;            % ori for punish
    results(trial_i,7) = Param.Stimuli.GratingOri(3);        % ori for control
 
    results(trial_i,8)  = results(trial_i,5) + (rand-0.5)*2*Param.Stimuli.OriJitter;
    results(trial_i,9)  = results(trial_i,6) + (rand-0.5)*2*Param.Stimuli.OriJitter;
    results(trial_i,10) = results(trial_i,7) + (rand-0.5)*2*Param.Stimuli.OriJitter;
     
    results(trial_i,11) = mod(trial_index(trial_i),3)+1; % 1 reward 2 punish 3 control
    results(trial_i,22) = results(trial_i,(results(trial_i,11)+1));  

    trial_onset = GetSecs;
    results(trial_i,13) = trial_onset;

    %% ITI
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
    vbl = Screen('Flip',wnd);
    results(trial_i,14) = vbl-trial_onset;

    %% stim
    [x,y]=meshgrid(-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize,-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize);
    angle1 = results(trial_i,8)/180*pi;   % reward
    angle2 = results(trial_i,9)/180*pi;   % punish
    angle3 = results(trial_i,10)/180*pi;  % control
    
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
    
    Screen('DrawTexture', wnd, mm1,[],[Param.Stimuli.Locations(results(trial_i,2),1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)+Param.Stimuli.OuterSize]);
    Screen('DrawTexture', wnd, mm2,[],[Param.Stimuli.Locations(results(trial_i,3),1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,3),2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,3),1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,3),2)+Param.Stimuli.OuterSize]);
    Screen('DrawTexture', wnd, mm3,[],[Param.Stimuli.Locations(results(trial_i,4),1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,4),2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,4),1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,4),2)+Param.Stimuli.OuterSize]);
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
    vbl = Screen('Flip',wnd,vbl+Param.Trial.ITI); %-Slack
    results(trial_i,15) = vbl-trial_onset-results(trial_i,14);
    
    %% delay
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
    vbl = Screen('Flip',wnd,vbl+Param.Trial.StiDura);     %-Slack
    results(trial_i,16) = vbl-trial_onset-sum(results(trial_i,14:15));
    
    %% report 
    [x,y]=meshgrid(-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize,-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize);
    phase = rand*2*pi;
    angle = rand*pi;
    results(trial_i,18) = angle*180/pi;
    
    sti4_final = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle)+x*cos(angle))+phase));
    sti4_final(sqrt(x.^2 + y.^2) > Param.Stimuli.OuterSize) = gray; %circle mask
    
    mm4 = Screen('MakeTexture', wnd, sti4_final);
    Screen('DrawTexture', wnd, mm4,[],[Param.Stimuli.Locations(results(trial_i,22),1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,22),2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,22),1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,22),2)+Param.Stimuli.OuterSize]);
    
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
    vbl = Screen('Flip',wnd,vbl+Param.Trial.Delay);%-Slack
    stimulus_onset = vbl;
    results(trial_i,17) = vbl-trial_onset-sum(results(trial_i,14:16));

    %% response
    is_true = 0;
    Angle_temp = 0;
    while (is_true == 0)
        [ifkey,RT_time,keyCode] = KbCheck;
        if keyCode(Param.Keys.Right)
            Angle_temp = Angle_temp+Param.Trial.Response;
        elseif keyCode(Param.Keys.Left)
            Angle_temp = Angle_temp-Param.Trial.Response;
        elseif keyCode(Param.Keys.Space)
            results(trial_i,21) = RT_time - stimulus_onset;  % RT
            is_true = 1;           
        elseif keyCode(Param.Keys.Esc)
            abort;
        end
    
        results(trial_i,12) = OriCircle180(Angle_temp+angle/pi*180);
    
            if results(trial_i,22) == results(trial_i,2)
               error = results(trial_i,12)-results(trial_i,8);

            elseif results(trial_i,22) == results(trial_i,3)
                error = results(trial_i,12)-results(trial_i,9);

            elseif results(trial_i,22) == results(trial_i,4)
                error = results(trial_i,12)-results(trial_i,10);
            end

            if error < (-90)
                results(trial_i,25) = error + 180;
            elseif error > 90
                results(trial_i,25) = error - 180;
            else
                results(trial_i,25) = error;
            end

            
        if ~is_true
            Screen('DrawTexture', wnd, mm4,[],[Param.Stimuli.Locations(results(trial_i,22),1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,22),2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,22),1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,22),2)+Param.Stimuli.OuterSize],Angle_temp);
            Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
            Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
            Screen('Flip',wnd);
        end
    end
    results(trial_i,20) = GetSecs - results(trial_i,13);
end

WaitSecs(2);

%% save the data
resultsDir = [CurrDir '\Results\test\' SubjID '\'];
if ~isdir(resultsDir)                                                       %isder和mkder可以用来创建新文件夹
mkdir(resultsDir);
end

cd(resultsDir);
results_name = [SubjID '_results_sess' num2str(sess_num) '_run' num2str(run_num) '.mat'];
save(results_name,'results','Param');
cd(CurrDir);
%% plot
hist(results(:,19));
S = std((results(:,19)));
M1 = mean(results(results(:,11)==1,19));
S1 = std((results(results(:,11)==1,19)));
M2 = mean(results(results(:,11)==2,19));
S2 = std((results(results(:,11)==2,19)));
M3 = mean(results(results(:,11)==3,19));
S3 = std((results(results(:,11)==3,19)));

%%
reset_test_gamma;
ShowCursor;
Screen('CloseAll');

delete *.asv
