function wm_motion_ver1(SubjID,session,block,MotionCoh)

parameter;

DrawFormattedText(wnd,'Press space to start','center','center', black);        
Screen('Flip',wnd);

is_true = 0;
while (is_true == 0)
    [~,~,keyCode] = KbCheck;
    if keyCode(Param.Keys.Space)
        is_true = 1;
    elseif keyCode(Param.Keys.EscPress)
        is_true = 1;
        abort;
    end
end

%***********Results************
%1：trial_i              %10：trial start time (fixation show time point)
%2：stim location        %11：生成运动刺激后时间点
%3：base direction       %12：生成运动刺激需要时间(<600ms)
%4：sample1 direction    %13：sample 1 show time
%5：sample2 direction    %14：ISI 1 show time
%6：cue                  %15：sample 2 show time
%7：subject response     %16：ISI 2 show time
%8：correct or not       %17：cue   show time
%9：responce time        %18：delay show time
%                        %19：test  呈现时间点
%                        %20：trial show time
%                        %21：coherence
%                        %22：test sample direction

results = zeros(Param.DisPrac.TrialNum,21);
trial_index = randperm(Param.DisPrac.TrialNum);
trial_index = mod(trial_index,Param.DisPrac.DirectionNum);

for trial_i = 1:Param.DisPrac.TrialNum
    
    %draw fixation
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);      
    vbl = Screen('Flip',wnd);
    results(trial_i,10) = vbl;
    
    results(trial_i,1) = trial_i;                             %1- No. trial
    Curr_cond = trial_index(trial_i)+1;   
    results(trial_i,2) = Param.DisPrac.StiLocation;           %2- not necessary if both sides
    results(trial_i,3) = Param.DisPrac.Directions(Curr_cond); %3- motion direction for stim1
    
    task_dif = Param.DisPrac.AngleDelta; 
    results(trial_i,21) = MotionCoh;
    
%     trial_order = ((rand-0.5)>=0);          % 0先呈现标准刺激 1后呈现标准刺激
    trial_delta_ori= sign(rand-0.5);        % -1: anti-clock 1: clock  
    if trial_delta_ori == 0
        trial_delta_ori = 1;
    end

    results(trial_i,4) = results(trial_i,3) + (rand-0.5)*2*Param.VarDir;
    results(trial_i,5) = results(trial_i,4) + trial_delta_ori * task_dif;

    MyDot_1 = Gen_DotMatrix_bm(results(trial_i,4),results(trial_i,21),Param.RDK.NumFrames,Param.RDK.DotNum,Param.RDK.InnerRadius,Param.RDK.OuterRadius,Param.RDK.StepPerMove);
    MyDot_2 = Gen_DotMatrix_bm(results(trial_i,5),results(trial_i,21),Param.RDK.NumFrames,Param.RDK.DotNum,Param.RDK.InnerRadius,Param.RDK.OuterRadius,Param.RDK.StepPerMove);    
    results(trial_i,11) = GetSecs; 
    results(trial_i,12) = results(trial_i,11)-results(trial_i,10); 

    % need add 'if results(trial_i,12)>600ms'?

    % display sample1
    for frame_i = 1:Param.RDK.NumFrames
        %%%motion
        Screen('DrawDots',wnd,MyDot_1{1,frame_i},Param.RDK.DotSize,black,Param.Stimuli.Locations(1,:),1);
        Screen('DrawDots',wnd,MyDot_1{1,frame_i},Param.RDK.DotSize,black,Param.Stimuli.Locations(results(trial_i,2),:),1);
        Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
        if frame_i == 1
            vbl = Screen('Flip',wnd,vbl+Param.Trial.ITI-Slack);
            results(trial_i,13)=vbl;  %%motion onset
        else
            vbl = Screen('Flip',wnd,vbl+Param.RDK.FramesPerMove/NominalFrameRate-Slack);
        end
    end
    
    
    % ISI1
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);    
    vbl = Screen('Flip',wnd, results(trial_i,13)+Param.RDK.Duration-Slack);
    results(trial_i,14)=vbl;
    
    
    % display sample2
    for frame_i = 1:Param.RDK.NumFrames
        %%%motion
        Screen('DrawDots',wnd,MyDot_2{1,frame_i},Param.RDK.DotSize,black,Param.Stimuli.Locations(1,:),1);
        Screen('DrawDots',wnd,MyDot_2{1,frame_i},Param.RDK.DotSize,black,Param.Stimuli.Locations(results(trial_i,2),:),1);
        Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
        if frame_i == 1
            vbl = Screen('Flip',wnd,vbl+Param.Trial.ISI-Slack);
            results(trial_i,15)=vbl;  %%motion onset
        else
            vbl = Screen('Flip',wnd,vbl+Param.RDK.FramesPerMove/NominalFrameRate-Slack);
        end
    end 
    
    % ISI2  
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);    
    vbl = Screen('Flip',wnd, results(trial_i,15)+Param.RDK.Duration-Slack);
    results(trial_i,16)=vbl;

    % set test
    sample_no = ((rand-0.5)>=0)+1;
    if sample_no ==1
        results(trial_i,22)=results(trial_i,4);
    else
        results(trial_i,22)=results(trial_i,5);
    end

    % display cue
    DrawFormattedText(wnd, Param.cue.Display(sample_no),'center','center', black);        
    vbl = Screen('Flip',wnd, results(trial_i,16)+Param.Trial.ISI-Slack);
    results(trial_i,17)=vbl;

    % delay
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
    vbl = Screen('Flip',wnd, results(trial_i,17)+Param.cue.Duration-Slack);
    results(trial_i,18)=vbl;

    % test sample
    MyDot_t = Gen_DotMatrix_bm(results(trial_i,22),results(trial_i,21),Param.RDK.TestNumFrames,Param.RDK.DotNum,Param.RDK.InnerRadius,Param.RDK.OuterRadius,Param.RDK.StepPerMove);
    for frame_i = 1:Param.RDK.TestNumFrames
        %%%motion
        Screen('DrawDots',wnd,MyDot_t{1,frame_i},Param.RDK.DotSize,black,Param.Stimuli.Locations(1,:),1);
        Screen('DrawDots',wnd,MyDot_t{1,frame_i},Param.RDK.DotSize,black,Param.Stimuli.Locations(results(trial_i,2),:),1);
        Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
        if frame_i == 1
            vbl = Screen('Flip',wnd,vbl+Param.delay.Duration-Slack);
            results(trial_i,19)=vbl;  %%motion onset
        else
            vbl = Screen('Flip',wnd,vbl+Param.RDK.FramesPerMove/NominalFrameRate-Slack);
        end
    end
    
    % response
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
    Screen('Flip',wnd, results(trial_i,19)+Param.RDK.TesDuration-Slack);    % i think there is no need to record this time point?
    is_true = 0;
    while (is_true == 0 && GetSecs - results(trial_i,19) < Param.Trial.MaxRT)
        [keyIsDown_1,RT_time,keyCode] = KbCheck;
        if keyCode(Param.Keys.Right) || keyCode(Param.Keys.two1) || keyCode(Param.Keys.two2)
            results(trial_i,7) = 2;              % subj response of whether test sample is same as sample2
            if sample_no == 2
                results(trial_i,8) = 1;          % response correct or not
            end
            results(trial_i,9) = RT_time - results(trial_i,19);  % 9反应时
            is_true = 1;
        elseif keyCode(Param.Keys.Left) || keyCode(Param.Keys.one1) || keyCode(Param.Keys.one2)
            results(trial_i,7) = 1;              % subj response of whether test sample is same as sample1
            if sample_no == 1
                results(trial_i,8) = 1;
            end
            results(trial_i,9) = RT_time - results(trial_i,19);
            is_true = 1;
        elseif keyCode(Param.Keys.EscPress)
            abort;
        end
        %没有反应的试次会被删除
    end
    
    if ~results(trial_i,8)                  
        beep;
    end
    
    %呈现刺激间隔  
    Screen('Drawtext',wnd,'');   
    Screen('Flip',wnd);
    while (GetSecs - results(trial_i,10) < Param.Trial.Duration)
    end
    results(trial_i,20) = GetSecs-results(trial_i,10);
end
    
prac_acc = sum(results(:,8))./Param.DisPrac.TrialNum;
disp('  ');
disp(['Accuracy: ' num2str(prac_acc)]);
disp('  ');

%数据存储
CurrDir = pwd;
resultsDir = [CurrDir '\Results\DisPrac\' SubjID '\'];
if ~isdir(resultsDir)                                                       %isder和mkder可以用来创建新文件夹
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
 
delete *.asv  