%% check parameters used
if exist([CurrDir '\Results\' SessName '\' SubjID '\' SubjID '_results_sess' num2str(Sess_Num) '_run' num2str(Run_Num) '.mat'],'file')
    disp(' ');
    disp([SubjID '_run' num2str(Run_Num) ' has been test, please enter a new run number!']); 
    disp(' ');
    abort;      
end

results = zeros(Param.Stairtest.StopRule1,15);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% results:
% 1 trial number              %11    ITI expectation
% 2 sti_loc                   %12    fix onset
% 3 sti_ori                   %13    ITI duration
% 4 angle of st1                 %14    st1 duration
% 5 angle of st2                 %15    ISI  
% 6 response                  %16    st2 duration    
% 7 task_diff                 
% 8 accuracy                   
% 9 reaction time             
%10 trial onset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% staircase settings
UD1 = PAL_AMUD_setupUD('up',Param.Stairtest.Up,'down',Param.Stairtest.Down);
UD1 = PAL_AMUD_setupUD(UD1,'StepSizeDown',Param.Stairtest.StepSizeDown,'StepSizeUp',Param.Stairtest.StepSizeUp,'stopcriterion',Param.Stairtest.StopCriterion,...
    'xMax',Param.Stairtest.xMax,'xMin',Param.Stairtest.xMin,'truncate','yes');

UD2 = PAL_AMUD_setupUD('up',Param.Stairtest.Up,'down',Param.Stairtest.Down);
UD2 = PAL_AMUD_setupUD(UD2,'StepSizeDown',Param.Stairtest.StepSizeDown,'StepSizeUp',Param.Stairtest.StepSizeUp,'stopcriterion',Param.Stairtest.StopCriterion,...
    'xMax',Param.Stairtest.xMax,'xMin',Param.Stairtest.xMin,'truncate','yes');

UD3 = PAL_AMUD_setupUD('up',Param.Stairtest.Up,'down',Param.Stairtest.Down);
UD3 = PAL_AMUD_setupUD(UD3,'StepSizeDown',Param.Stairtest.StepSizeDown,'StepSizeUp',Param.Stairtest.StepSizeUp,'stopcriterion',Param.Stairtest.StopCriterion,...
    'xMax',Param.Stairtest.xMax,'xMin',Param.Stairtest.xMin,'truncate','yes');

if Run_Num == 1 | Run_Num == 2
    UD1 = PAL_AMUD_setupUD(UD1,'startvalue', Param.Stairtest.start,'stoprule',Param.Stairtest.StopRule1);
    UD2 = PAL_AMUD_setupUD(UD2,'startvalue', Param.Stairtest.start,'stoprule',Param.Stairtest.StopRule1);
    UD3 = PAL_AMUD_setupUD(UD3,'startvalue', Param.Stairtest.start,'stoprule',Param.Stairtest.StopRule1);
else 
    previou_results = load([CurrDir '\Results\' SessName '\' SubjID '\' SubjID '_results_sess' num2str(Sess_Num)  '_run' num2str(Run_Num-1) '.mat']);
    UD1 = PAL_AMUD_setupUD(UD1,'startvalue', Param.Stairtest.start,'stoprule',Param.Stairtest.StopRule2);
    UD2 = PAL_AMUD_setupUD(UD2,'startvalue', Param.Stairtest.start,'stoprule',Param.Stairtest.StopRule2);
    UD3 = PAL_AMUD_setupUD(UD3,'startvalue', Param.Stairtest.start,'stoprule',Param.Stairtest.StopRule2);
end

%% start
DrawFormattedText(wnd,'Press space to start!','center','center', black);
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
pause(0.5);


%% Go!
while (~UD1.stop) | (~UD2.stop) | (~UD3.stop)
    trial_i = size(UD1.response,2)+size(UD2.response,2)+size(UD3.response,2)+1;
    
    UD_continue = [];
    if (~UD1.stop)
        UD_continue = [UD_continue,1];
    end

    if (~UD2.stop)
        UD_continue = [UD_continue,2];
    end

    if (~UD3.stop)
        UD_continue = [UD_continue,3];
    end

    results(trial_i,1) = trial_i;
    results(trial_i,2) = location_used;

    if Sess_Num == 1
        ori_temp = Shuffle(UD_continue);
        ori_temp = ori_temp(1);
    elseif Sess_Num == 2
        ori_temp = mod((trial_i-1),3)+1;
    end

    results(trial_i,3) = Param.Stimuli.GratingOri(ori_temp);

    if ori_temp == 1
        results(trial_i,7) = UD1.xCurrent;
    elseif ori_temp == 2
        results(trial_i,7) = UD2.xCurrent;
    else 
        results(trial_i,7) = UD3.xCurrent;
    end
           

    trial_onset = GetSecs;
    results(trial_i,10) = trial_onset; 

    trial_delta_ori = ((rand-0.5)>=0)*2-1;  
    trial_order = ((rand-0.5)>=0);          
    if trial_order    %1ref first %0 test first
        results(trial_i,4) = results(trial_i,3) + (rand-0.5)*2*Param.Stimuli.OriJitter;
        results(trial_i,5) = results(trial_i,4) + trial_delta_ori*results(trial_i,7);
    else
        results(trial_i,5) = results(trial_i,3) + (rand-0.5)*2*Param.Stimuli.OriJitter;
        results(trial_i,4) = results(trial_i,5) + trial_delta_ori*results(trial_i,7);    
    end

    %% ITI
    rhyme = randperm(length(Param.Trial.ITI),1);
    results(trial_i,11) = Param.Trial.ITI(rhyme);

    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
    vbl = Screen('Flip',wnd);
    results(trial_i,12) = vbl-trial_onset; 
    
    %% sti1
    [x,y]=meshgrid(-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize,-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize);
    phase = rand*2*pi; 
    angle = results(trial_i,4)/180*pi;
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

    Screen('DrawTexture', wnd, mm,[],[Param.Stimuli.Locations(results(trial_i,2),1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)+Param.Stimuli.OuterSize]);
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
    vbl = Screen('Flip',wnd,vbl+Param.Trial.ITI(rhyme));
    results(trial_i,13) = vbl-trial_onset-results(trial_i,12);

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

    % sin grating
    Gabor_sti = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle)+x*cos(angle))+phase).*exp(-(x.^2 + y.^2)/(2.*Param.Stimuli.SmoothSD.^2))); 
    Gabor_sti(sqrt(x.^2 + y.^2) > Param.Stimuli.OuterSize) = gray;
    
    Gabor_sti_final = repmat(Gabor_sti,[1,1,3]); 
    Gabor_sti_final(:,:,4) = annulus*white;
    
    mm = Screen('MakeTexture', wnd, Gabor_sti_final); 
    Screen('DrawTexture', wnd, mm,[],[Param.Stimuli.Locations(results(trial_i,2),1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)+Param.Stimuli.OuterSize]);
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
    vbl = Screen('Flip',wnd,vbl+Param.Trial.ISI); %-Slack
    stimulus_onset = vbl;
    results(trial_i,15) = vbl-trial_onset-sum(results(trial_i,12:14));

    %% response collection
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
    vbl = Screen('Flip',wnd,vbl+Param.Trial.StiDura);%-Slack
    results(trial_i,16) = vbl-trial_onset-sum(results(trial_i,12:15));      
    
    is_true = 0;
    while (is_true == 0 & GetSecs - stimulus_onset < Param.Trial.MaxRT)
        [ifkey,RT_time,keyCode] = KbCheck;
        if keyCode(Param.Keys.Right)
            results(trial_i,6) = 2;              %6 sti in interval 2
            if results(trial_i,5) > results(trial_i,4)
                results(trial_i,8) = 1;          %8 correct or not
            end
            results(trial_i,9) = RT_time - stimulus_onset;  %9 RT
            is_true = 1;
        elseif keyCode(Param.Keys.Left)
            results(trial_i,6) = 1;              
            if results(trial_i,5) < results(trial_i,4)
                results(trial_i,8) = 1;
            end
            results(trial_i,9) = RT_time - stimulus_onset;
            is_true = 1;
        elseif keyCode(Param.Keys.Esc)
            abort;
        end
    end

    if ~results(trial_i,8)
%                     beep;
        DrawFormattedText(wnd,'Wrong!','center','center', black);        
    else
        DrawFormattedText(wnd,'Correct!','center','center', black);        
    end
    Screen('Flip',wnd);


    if ori_temp == 1
        UD1 = PAL_AMUD_updateUD(UD1, results(trial_i,8)); %update UD structure
    elseif ori_temp == 2
        UD2 = PAL_AMUD_updateUD(UD2, results(trial_i,8)); %update UD structure
    else
        UD3 = PAL_AMUD_updateUD(UD3, results(trial_i,8)); %update UD structure
    end
        
    while (GetSecs - trial_onset < Param.Trial.Duration)
    end         
end
 
% take a break (awake consolidation)

       if (Run_Num > 3) && (mod(Run_Num,2)==0)
           Screen('FillRect',wnd,black,Param.Settings.ScrnResolution);
           Screen('Flip',wnd);
           pause(240);
           Screen('FillRect',wnd,gray,Param.Settings.ScrnResolution);
           sample_ratio = 0;
           sample;
           pause(10);
       end

%% save the data

threshold_value1 = PAL_AMUD_analyzeUD(UD1, 'trials', Param.Stairtest.trialthresh);
threshold_value2 = PAL_AMUD_analyzeUD(UD2, 'trials', Param.Stairtest.trialthresh);
threshold_value3 = PAL_AMUD_analyzeUD(UD3, 'trials', Param.Stairtest.trialthresh);

resultsDir = [CurrDir '\Results\' SessName '\' SubjID '\'];
if ~isdir(resultsDir)                                                       
    mkdir(resultsDir);
end

cd(resultsDir);
results_name = [SubjID '_results_sess' num2str(Sess_Num) '_run' num2str(Run_Num) '.mat'];
save(results_name,'results','UD1','UD2','UD3','threshold_value1','threshold_value2','threshold_value3','Param');
cd(CurrDir); 

%% plot
for figure_i = 1:length(Param.Stimuli.GratingOri)
    figure(figure_i);
    if figure_i == 1
        end_trial = size(UD1.x,2);
        task_diff_temp = UD1.x;
    elseif figure_i == 2
        end_trial = size(UD2.x,2);
        task_diff_temp = UD2.x;
    else
        end_trial = size(UD3.x,2);
        task_diff_temp = UD3.x;
    end
   
    plot(1:end_trial,task_diff_temp(1:end_trial));
    axis([0 Param.Stairtest.MaxTrial 0 15]);  
end
%%
% reset_test_gamma;
ShowCursor;
Screen('CloseAll');
 
delete *.asv