% orientation detection task
% quest

%% check parameters used
if exist([CurrDir '\Results\' SessName '\' SubjID '\' SubjID '_results_sess' num2str(Sess_Num)  '_run' num2str(Run_Num) '.mat'],'file')
    disp(' ');
    disp([SubjID  '_run' num2str(Run_Num) ' has been test, please enter a new run number!']); 
    disp(' ');
    abort;      
end

results = zeros(Param.RF.stoprule*length(Param.Stimuli.GratingOri),18);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% results:
% 1 trial number              %11    ITI expectation
% 2 sti_loc                   %12    fix onset
% 3 sti_ori                   %13    ITI duration
% 4 SN of st1                 %14    cue duration
% 5 SN of st2                 %15    gap duration  
% 6 response                  %16    st1 duration    
% 7 task_diff                 %17    ISI
% 8 accuracy                  %18    st2 duration 
% 9 reaction time             %19    tone
%10 trial onset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% quest settings
RF1 = PAL_AMRF_setupRF('priorAlphaRange', Param.RF.alphas,...
    'stopcriterion',Param.RF.stopcriterion,'beta',Param.RF.beta,...
    'gamma',Param.RF.gamma,'lambda',Param.RF.lambda,'PF',Param.RF.PFfit,'meanmode',Param.RF.meanmode);
RF1 = PAL_AMRF_setupRF(RF1,'xMax',Param.RF.xMax,'xMin',Param.RF.xMin);

RF2 = PAL_AMRF_setupRF('priorAlphaRange', Param.RF.alphas,...
    'stopcriterion',Param.RF.stopcriterion,'beta',Param.RF.beta,...
    'gamma',Param.RF.gamma,'lambda',Param.RF.lambda,'PF',Param.RF.PFfit,'meanmode',Param.RF.meanmode);
RF2 = PAL_AMRF_setupRF(RF2,'xMax',Param.RF.xMax,'xMin',Param.RF.xMin);

RF3 = PAL_AMRF_setupRF('priorAlphaRange', Param.RF.alphas,...
    'stopcriterion',Param.RF.stopcriterion,'beta',Param.RF.beta,...
    'gamma',Param.RF.gamma,'lambda',Param.RF.lambda,'PF',Param.RF.PFfit,'meanmode',Param.RF.meanmode);
RF3 = PAL_AMRF_setupRF(RF3,'xMax',Param.RF.xMax,'xMin',Param.RF.xMin);

RF1 = PAL_AMRF_setupRF(RF1,'prior', Param.RF.prior,'stoprule',Param.RF.stoprule);
RF2 = PAL_AMRF_setupRF(RF2,'prior', Param.RF.prior,'stoprule',Param.RF.stoprule);
RF3 = PAL_AMRF_setupRF(RF3,'prior', Param.RF.prior,'stoprule',Param.RF.stoprule);

%start   
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
pause(0.2);

%% Go!
while (~RF1.stop) || (~RF2.stop) || (~RF3.stop)
    trial_i = size(RF1.response,2)+size(RF2.response,2)+size(RF3.response,2)+1;

    RF_continue = [];
    if (~RF1.stop)
        RF_continue = [RF_continue,1];
    end
    
    if (~RF2.stop)
        RF_continue = [RF_continue,2];
    end
    
    if (~RF3.stop)
        RF_continue = [RF_continue,3];
    end
    results(trial_i,1) = trial_i;
    results(trial_i,2) = location_used;

    if Sess_Num == 1
        ori_temp = Shuffle(RF_continue);
        ori_temp = ori_temp(1);
    elseif Sess_Num == 2
        ori_temp = mod((trial_i-1),3)+1;
    end
    results(trial_i,3) = Param.Stimuli.GratingOri(ori_temp);
  
    if ori_temp == 1
        results(trial_i,7) = RF1.xCurrent;
    elseif ori_temp == 2
        results(trial_i,7) = RF2.xCurrent;
    else
        results(trial_i,7) = RF3.xCurrent;
    end
           
    trial_onset = GetSecs;
    results(trial_i,10) = trial_onset; 

    %% ITI
    rhyme = randperm(length(Param.Trial.ITI),1);
    results(trial_i,11) = Param.Trial.ITI(rhyme);

    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
    vbl = Screen('Flip',wnd);
    results(trial_i,12) = vbl-trial_onset; 
    
    %% gen stimuli
    % noise pattern 
    [x,y]=meshgrid(-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize,-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize);
    phase = rand*2*pi; 
    angle = rand*pi;
    Noise_sti1 = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle)+x*cos(angle))+phase));
    
    Noise_sti1 = Matrix_shuffle(Noise_sti1);
    Noise_sti1(sqrt(x.^2 + y.^2) > Param.Stimuli.OuterSize) = gray; 
    
    % signal pattern
    [x,y]=meshgrid(-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize,-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize);
    phase = rand*2*pi; 
    angle = rand*pi;
    Noise_sti2 = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle)+x*cos(angle))+phase)); 
    Noise_sti2 = Matrix_shuffle(Noise_sti2);
    Noise_sti2(sqrt(x.^2 + y.^2) > Param.Stimuli.OuterSize) = gray; 

    [x,y]=meshgrid(-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize,-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize);
    phase = rand*2*pi; 
    angle = results(trial_i,3)/180*pi;
    Signal_sti = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle)+x*cos(angle))+phase).*exp(-(x.^2 + y.^2)/2/(Param.Stimuli.SmoothSD.^2))); 
    Signal_sti(sqrt(x.^2 + y.^2) > Param.Stimuli.OuterSize) = gray;
   
    % selection
    overlay_temp = ones(size(x));
    overlay_temp(sqrt(x.^2 + y.^2) > Param.Stimuli.OuterSize) = 0;
    overlay_temp(sqrt(x.^2 + y.^2) < Param.Stimuli.InnerSize) = 0;
    selection_temp = find(overlay_temp==1);
    overlay_temp = zeros(size(x));

    signal_pixels_num = round(results(trial_i,7) * length(selection_temp));
    rand_temp = randperm(length(selection_temp));   
    signal_pixels = rand_temp(1:signal_pixels_num);
    overlay_temp(selection_temp(signal_pixels)) = 1;

    % add together 
    trial_order = ((rand-0.5)>=0);          % 0signal first 1noise first
    if trial_order
        results(trial_i,5) = results(trial_i,7);
        results(trial_i,4) = 0;
        Sti1_final = Noise_sti1;
        Sti2_final = Noise_sti2.*(1-overlay_temp)+Signal_sti.*overlay_temp;
    else
        results(trial_i,4) = results(trial_i,7);
        results(trial_i,5) = 0;
        Sti1_final = Noise_sti2.*(1-overlay_temp)+Signal_sti.*overlay_temp;
        Sti2_final = Noise_sti1;         
    end
    
    %% Cue   
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
    vbl = Screen('Flip',wnd,vbl+Param.Trial.ITI(rhyme));
    results(trial_i,13) = vbl-trial_onset-results(trial_i,12);
    
    %Tag
    if Sess_Num == 1
        Snd ('Play',beep0);
% PsychPortAudio('FillBuffer', pahandle, [beep0;beep0]);
% PsychPortAudio('Start', pahandle, 1,0,1);

    elseif Sess_Num == 2
        Snd ('Play', eval(['beep' num2str(ori_temp)]));
%         Beeper(frequency, [fVolume], [durationSec]);
    end
    results(trial_i,19) = ori_temp;
    
    %% Gap   
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
    vbl = Screen('Flip',wnd,vbl+Param.Trial.Cue);
    results(trial_i,14) = vbl-trial_onset-sum(results(trial_i,12:13));
    
    %% sti1
    mm = Screen('MakeTexture', wnd, Sti1_final);
    Screen('DrawTexture', wnd, mm,[],[Param.Stimuli.Locations(results(trial_i,2),1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)+Param.Stimuli.OuterSize]);
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
    vbl = Screen('Flip',wnd,vbl+Param.Trial.Gap); %-Slack
    results(trial_i,15) = vbl-trial_onset-sum(results(trial_i,12:14));  
    
    %% ISI
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
    vbl = Screen('Flip',wnd,vbl+Param.Trial.StiDura); %-Slack
    results(trial_i,16) = vbl-trial_onset-sum(results(trial_i,12:15));
    
    %% sti2
    mm = Screen('MakeTexture', wnd, Sti2_final); 
    Screen('DrawTexture', wnd, mm,[],[Param.Stimuli.Locations(results(trial_i,2),1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)+Param.Stimuli.OuterSize]);
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
    vbl = Screen('Flip',wnd,vbl+Param.Trial.ISI); %-Slack
    stimulus_onset = vbl;
    results(trial_i,17) = vbl-trial_onset-sum(results(trial_i,12:16));
 
    %% response collection
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
    vbl = Screen('Flip',wnd,vbl+Param.Trial.StiDura);%-Slack
    results(trial_i,18) = vbl-trial_onset-sum(results(trial_i,12:17));      
    
    is_true = 0;
    while (is_true == 0 & GetSecs - stimulus_onset < Param.Trial.MaxRT)
        [ifkey,RT_time,keyCode] = KbCheck;
        if keyCode(Param.Keys.J)
            results(trial_i,6) = 2;              %6 clockwise
            if results(trial_i,5)> results(trial_i,4)
                results(trial_i,8) = 1;          %8 correct or not
            end
            results(trial_i,9) = RT_time - stimulus_onset;  %9 RT
            is_true = 1;
        elseif keyCode(Param.Keys.F)
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
        %             beep;
        DrawFormattedText(wnd,'Wrong!','center','center', black);        
    else
        DrawFormattedText(wnd,'Correct!','center','center', black);        
    end
    Screen('Flip',wnd);

    if ori_temp == 1
        RF1 = PAL_AMRF_updateRF(RF1,results(trial_i,7), results(trial_i,8));
    elseif ori_temp == 2
        RF2 = PAL_AMRF_updateRF(RF2,results(trial_i,7), results(trial_i,8));
    elseif ori_temp == 3
        RF3 = PAL_AMRF_updateRF(RF3,results(trial_i,7), results(trial_i,8));
    end

    while (GetSecs - trial_onset < Param.Trial.Duration)
    end
end

%% save the data
resultsDir = [CurrDir '\Results\' SessName '\' SubjID '\'];
if ~isdir(resultsDir)                                                   
    mkdir(resultsDir);
end

cd(resultsDir);
curr_threshold1 = RF1.mean;
curr_threshold2 = RF2.mean;
curr_threshold3 = RF3.mean;

results_name = [SubjID '_results_sess' num2str(Sess_Num) '_run' num2str(Run_Num) '.mat'];
save(results_name,'results','RF1','RF2','RF3','Param','curr_threshold1','curr_threshold2','curr_threshold3');
cd(CurrDir); 

%% plot 
for figure_i = 1:length(Param.Stimuli.GratingOri)
    figure(figure_i); 
    if figure_i == 1
        trial_num_temp = 1:length(RF1.x);    
        plot(trial_num_temp,RF1.x,'k');
        hold on;
        plot(trial_num_temp(RF1.response == 1),RF1.x(RF1.response == 1),'ko', 'MarkerFaceColor','k');
        plot(trial_num_temp(RF1.response == 0),RF1.x(RF1.response == 0),'ko', 'MarkerFaceColor','w');       
        axis([0 max(trial_num_temp)+1 0 max(RF1.x)+(max(RF1.x)-min(RF1.x))/10]);     
    elseif figure_i == 2
        trial_num_temp = 1:length(RF2.x); 
        plot(trial_num_temp,RF2.x,'k');
        hold on;
        plot(trial_num_temp(RF2.response == 1),RF2.x(RF2.response == 1),'ko', 'MarkerFaceColor','k');
        plot(trial_num_temp(RF2.response == 0),RF2.x(RF2.response == 0),'ko', 'MarkerFaceColor','w');   
        axis([0 max(trial_num_temp)+1 0 max(RF2.x)+(max(RF2.x)-min(RF2.x))/10]);  
    elseif figure_i == 3
        trial_num_temp = 1:length(RF3.x); 
        plot(trial_num_temp,RF3.x,'k');
        hold on;
        plot(trial_num_temp(RF3.response == 1),RF3.x(RF3.response == 1),'ko', 'MarkerFaceColor','k');
        plot(trial_num_temp(RF3.response == 0),RF3.x(RF3.response == 0),'ko', 'MarkerFaceColor','w');       
        axis([0 max(trial_num_temp)+1 0 max(RF3.x)+(max(RF3.x)-min(RF3.x))/10]);       
    end 

    set(gca,'FontSize',16);
    xlabel('Trials');
    ylabel('Difference');  
end

reset_test_gamma;
ShowCursor;
Screen('CloseAll');
 
delete *.asv