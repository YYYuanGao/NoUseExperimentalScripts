function wm_motion_ver3(SubjID,session,block,MotionCoh)

parameter;

DrawFormattedText(wnd,'Press space to start','center','center', black);        
Screen('Flip',wnd);

is_true = 0;
while (is_true == 0)
    [~,~,keyCode] = KbCheck;
    if keyCode(Param.Keys.Space)
        is_true = 1;
    elseif keyCode(Param.Keys.EscPress)
        abort;
        error('Aborting due to user cancellation...');
    end
end

%*************************************
%***********Results Matrix************
%*************************************
% 1：trial_i                 %11：trial(fixation bg) start time point
% 2：motion coherence        %12：time point generated motion dots
% 3：base direction          %13：duration to generate motion dots <600ms
% 4：sample1 direction       %14：time point to display sample1
% 5：sample2 direction       %15：time point to display ISI1
% 6：cue no.                 %16：time point to display sample1
% 7：test no.                %17：time point to display ISI2
% 8：test direction          %18：time point to display cue
% 9：subject's response      %19：time point to display delay
% 10:correct or not         %20：time point to display test sample
%                           %21：reaction time
%                           %22：trial duration

results = zeros(Param.DisPrac.TrialNum,21);
trial_index = randperm(Param.DisPrac.TrialNum);
trial_index = mod(trial_index,Param.DisPrac.DirectionNum);
if mod(Param.DisPrac.TrialNum,Param.DisPrac.DirectionNum)~=0
    trial_index = randperm(Param.DisPrac.TrialNum-mod(Param.DisPrac.TrialNum,Param.DisPrac.DirectionNum));
    trial_index = mod(trial_index,Param.DisPrac.DirectionNum);
    mod_trial = [1,0,3];
    if mod(block,2)
        mod_trial = [4,5,2];
    end
    trial_index = cat(2,trial_index,mod_trial);
end
%% 
for trial_i = 1:Param.DisPrac.TrialNum
    
    %draw fixation
    Screen('FillOval', wnd, Param.Fixation.OvalColor, Param.Fixation.OvalLoc);
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);      
    vbl = Screen('Flip',wnd);
    results(trial_i,11) = vbl;
    
    results(trial_i,1) = trial_i;                               %1- No. trial
    Curr_cond = trial_index(trial_i)+1;
    stimLoc = Param.DisPrac.StiLocation;
    results(trial_i,2) = MotionCoh;                             %2- motion coherence
    results(trial_i,3) = Param.DisPrac.Directions(Curr_cond);   %3- baseline motion direction

    % set sample directions, baseline and jitter+baseline
    smp_order = ((rand-0.5)>=0);
    if smp_order
        results(trial_i,4) = results(trial_i,3) + (rand-0.5)*2*Param.VarDir;    % AngleDerta or VarDir
%       results(trial_i,5) = results(trial_i,4) + trial_delta_ori * task_dif;
        results(trial_i,5) = results(trial_i,3);
    else
        results(trial_i,4) = results(trial_i,3);    % AngleDerta or VarDir
        results(trial_i,5) = results(trial_i,3) + (rand-0.5)*2*Param.VarDir;
    end

    MyDot_1 = Gen_DotMatrix_bm(results(trial_i,4),results(trial_i,2),Param.RDK.NumFrames,Param.RDK.DotNum,Param.RDK.InnerRadius,Param.RDK.OuterRadius,Param.RDK.StepPerMove);
    MyDot_2 = Gen_DotMatrix_bm(results(trial_i,5),results(trial_i,2),Param.RDK.NumFrames,Param.RDK.DotNum,Param.RDK.InnerRadius,Param.RDK.OuterRadius,Param.RDK.StepPerMove);    
    results(trial_i,12) = GetSecs; 
    results(trial_i,13) = results(trial_i,12)-results(trial_i,11); 

    % need add 'if results(trial_i,12)>600ms'?

    % display sample1
    for frame_i = 1:Param.RDK.NumFrames
        %%%motion
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
    
    % ISI1
    Screen('FillOval', wnd, Param.Fixation.OvalColor, Param.Fixation.OvalLoc);
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);    
    vbl = Screen('Flip',wnd, results(trial_i,14)+Param.RDK.Duration-Slack);
    results(trial_i,15)=vbl;
    
    
    % display sample2
    for frame_i = 1:Param.RDK.NumFrames
        %%%motion
%         Screen('DrawDots',wnd,MyDot_2{1,frame_i},Param.RDK.DotSize,black,Param.Stimuli.Locations(1,:),1);
        Screen('DrawDots',wnd,MyDot_2{1,frame_i},Param.RDK.DotSize,black,Param.Stimuli.Locations(stimLoc,:),1);
        Screen('FillOval', wnd, Param.Fixation.OvalColor, Param.Fixation.OvalLoc);
        Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
        if frame_i == 1
            vbl = Screen('Flip',wnd,vbl+Param.Trial.ISI-Slack);
            results(trial_i,16)=vbl;  %%motion onset
        else
            vbl = Screen('Flip',wnd,vbl+Param.RDK.FramesPerMove/NominalFrameRate-Slack);
        end
    end 
    
    % ISI2  
    Screen('FillOval', wnd, Param.Fixation.OvalColor, Param.Fixation.OvalLoc);
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);    
    vbl = Screen('Flip',wnd, results(trial_i,16)+Param.RDK.Duration-Slack);
    results(trial_i,17)=vbl;

    % set test
    cue_no = ((rand-0.5)>=0)+1;
    results(trial_i,6)=cue_no;

    % display cue
    DrawFormattedText(wnd, Param.cue.Display(cue_no),'center','center', black);        
    vbl = Screen('Flip',wnd, results(trial_i,17)+Param.Trial.ISI-Slack);
    results(trial_i,18)=vbl;

    % delay   
    Screen('FillOval', wnd, Param.Fixation.OvalColor, Param.Fixation.OvalLoc);
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
    vbl = Screen('Flip',wnd, results(trial_i,18)+Param.cue.Duration-Slack);
    results(trial_i,19)=vbl;

    % set test
    test_no = ((rand-0.5)>=0)+1;
    results(trial_i,7) = test_no;
    if test_no ==1
        results(trial_i,8)=results(trial_i,4);
    else
        results(trial_i,8)=results(trial_i,5);
    end

    % test sample
    MyDot_t = Gen_DotMatrix_bm(results(trial_i,8),results(trial_i,2),Param.RDK.TestNumFrames,Param.RDK.DotNum,Param.RDK.InnerRadius,Param.RDK.OuterRadius,Param.RDK.StepPerMove);
    for frame_i = 1:Param.RDK.TestNumFrames
        %%%motion
%         Screen('DrawDots',wnd,MyDot_t{1,frame_i},Param.RDK.DotSize,black,Param.Stimuli.Locations(1,:),1);
        Screen('DrawDots',wnd,MyDot_t{1,frame_i},Param.RDK.DotSize,black,Param.Stimuli.Locations(stimLoc,:),1);
        Screen('FillOval', wnd, Param.Fixation.OvalColor, Param.Fixation.OvalLoc);
        Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
        if frame_i == 1
            vbl = Screen('Flip',wnd,vbl+Param.delay.Duration-Slack);
            results(trial_i,20)=vbl;  %%motion onset
        else
            vbl = Screen('Flip',wnd,vbl+Param.RDK.FramesPerMove/NominalFrameRate-Slack);
        end
    end
    
    % response
    Screen('FillOval', wnd, Param.Fixation.OvalColor, Param.Fixation.OvalLoc);
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
    Screen('Flip',wnd, results(trial_i,20)+Param.RDK.TesDuration-Slack);    % i think there is no need to record this time point?
    is_true = 0;

    brk_flag = 0;
    while (is_true == 0 && GetSecs - results(trial_i,20) < Param.Trial.MaxRT)
        [keyIsDown_1,RT_time,keyCode] = KbCheck;
        if keyCode(Param.Keys.Right) || keyCode(Param.Keys.two1) || keyCode(Param.Keys.two2)
            results(trial_i,9) = 2;              % subj responds the test sample is not corresponding to cue
            if test_no ~= cue_no
                results(trial_i,10) = 1;          % Correct: test is different to cue
            end
            results(trial_i,21) = RT_time - results(trial_i,20);  % 21- reaction time
            is_true = 1;
        elseif keyCode(Param.Keys.Left) || keyCode(Param.Keys.one1) || keyCode(Param.Keys.one2)
            results(trial_i,9) = 1;              % subj responds the test sample is corresponding to cue
            if test_no == cue_no
                results(trial_i,10) = 1;
            end
            results(trial_i,21) = RT_time - results(trial_i,20);
            is_true = 1;
        elseif keyCode(Param.Keys.EscPress)
            abort;
            brk_flag = 1;
            break
        end
        %trials without response will be deleted
    end
    if brk_flag
        break
    end
    
    if ~results(trial_i,10)                  
        beep;
    end
    
    %ITI  
%     Screen('Drawtext',wnd,'');   
    Screen('FillOval', wnd, Param.Fixation.OvalColor, Param.Fixation.OvalLoc);  
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
    Screen('Flip',wnd);
    while (GetSecs - results(trial_i,11) < Param.Trial.Duration-0.5)
    end
    Screen('FillOval', wnd, Param.Fixation.CrossColor, Param.Fixation.OvalLoc);  
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.OvalColor, [], 1);
    Screen('Flip',wnd);
    while (GetSecs - results(trial_i,11) < Param.Trial.Duration)
    end
    results(trial_i,22) = GetSecs-results(trial_i,11);
end

%%
prac_acc = sum(results(:,10))./Param.DisPrac.TrialNum;
disp('  ');
disp(['Accuracy: ' num2str(prac_acc)]);
disp('  ');

% save data
CurrDir = pwd;
resultsDir = [CurrDir '\Results\DisPrac\' SubjID '\'];
if ~isdir(resultsDir)
    mkdir(resultsDir);
end

cd(resultsDir);
results_name = [SubjID '_DisPrac_results_session' num2str(session) '_block' num2str(block) '.mat'];
save(results_name,'results','prac_acc');
cd(CurrDir); 

warning on;
%reset_test_gamma;
ShowCursor;
Screen('CloseAll');

if brk_flag
    error('Aborting due to user cancellation...');
end
 
delete *.asv  

