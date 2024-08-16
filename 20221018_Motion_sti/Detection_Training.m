function Detection_Training(SubjID,session,block,MotionDirection,MotionCoh)

% subj_id = 'jiake'
% by Ke Jia
% Last modified 2021/10/28 10:15

parameter;

%*************************************
%*************Exp Start***************
%*************************************
DrawFormattedText(wnd,'Press space to start','center','center', black);        
Screen('Flip',wnd);

is_true = 0;
while (is_true == 0)
    [~,~,keyCode] = KbCheck;
    if keyCode(Param.Keys.Space)
        is_true = 1;
    elseif keyCode(Param.Keys.EscPress)
        abort;
        is_true = 1;
    end
end


%*************************************
%***********Results Matrix************
%*************************************
%1：trial_i       %10：试次开始时间=fixation呈现时间点
%2：刺激位置       %11：生成运动刺激后时间点
%3：运动方向       %12：生成运动刺激需要时间,注意这个值要小于600ms
%4：一致性一       %13：运动刺激1呈现时间点
%5：一致性二       %14：刺激间隔呈现时间点
%6：被试反应       %15：运动刺激2呈现时间点
%7：一致水平       %16：fixation呈现时间点
%8：正确与否       %17：试次长度
%9：反应时间

results = zeros(Param.Staircase.MaxTrial,17);

task_dif = MotionCoh;
UD = PAL_AMUD_setupUD('up',Param.Staircase.Up,'down',Param.Staircase.Down);
UD = PAL_AMUD_setupUD(UD,'StepSizeDown',Param.Staircase.StepSizeDownDetect,'StepSizeUp',Param.Staircase.StepSizeUpDetect,'stopcriterion',Param.Staircase.StopCriterion,'stoprule',Param.Staircase.StopRule,'startvalue',task_dif,...
        'xMax',Param.Staircase.xMaxDetect,'xMin',Param.Staircase.xMin,'truncate','yes');

while (~UD.stop & (size(UD.response)<Param.Staircase.MaxTrial)) 
    trial_i = size(UD.response,2)+1;
    results(trial_i,1) = trial_i;                                %1试次数
    
    %draw fixation
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);      
    vbl = Screen('Flip',wnd);
    results(trial_i,10) = vbl;
      
    results(trial_i,2) = Param.DetectTrain.StiLocation;           %2刺激位置
    results(trial_i,3) = MotionDirection;                        %3运动方向
    results(trial_i,7) = UD.xCurrent;
    
    trial_order = ((rand-0.5)>=0);          % 0先呈现标准刺激 1后呈现标准刺激
    if trial_order
        results(trial_i,5) = results(trial_i,7);  
        results(trial_i,4) = 0;    
    else
        results(trial_i,4) = results(trial_i,7); 
        results(trial_i,5) = 0;            
    end

    MyDot_1 = Gen_DotMatrix_bm(results(trial_i,3),results(trial_i,4),Param.RDK.NumFrames,Param.RDK.DotNum,Param.RDK.InnerRadius,Param.RDK.OuterRadius,Param.RDK.StepPerMove);
    MyDot_2 = Gen_DotMatrix_bm(results(trial_i,3),results(trial_i,5),Param.RDK.NumFrames,Param.RDK.DotNum,Param.RDK.InnerRadius,Param.RDK.OuterRadius,Param.RDK.StepPerMove);    
    results(trial_i,11) = GetSecs; 
    results(trial_i,12) = results(trial_i,11)-results(trial_i,10); 

    %呈现运动点刺激
    for frame_i = 1:Param.RDK.NumFrames
        %%%motion
        Screen('DrawDots',wnd,MyDot_1{1,frame_i},Param.RDK.DotSize,black,Param.Stimuli.Locations(results(trial_i,2),:),1);
        Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
        if frame_i == 1
            vbl = Screen('Flip',wnd,vbl+Param.Trial.ITI-Slack);
            results(trial_i,13)=vbl;  %%motion onset
        else
            vbl = Screen('Flip',wnd,vbl+Param.RDK.FramesPerMove/NominalFrameRate-Slack);
        end
    end 
    
    
    %呈现刺激间隔  
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);    
    vbl = Screen('Flip',wnd, results(trial_i,13)+Param.RDK.Duration-Slack);
    results(trial_i,14)=vbl;
    
    
    %呈现运动点刺激
    for frame_i = 1:Param.RDK.NumFrames
        %%%motion
        Screen('DrawDots',wnd,MyDot_2{1,frame_i},Param.RDK.DotSize,black,Param.Stimuli.Locations(results(trial_i,2),:),1);
        Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
        if frame_i == 1
            vbl = Screen('Flip',wnd,vbl+Param.Trial.ISI-Slack);
            results(trial_i,15)=vbl;  %%motion onset
        else
            vbl = Screen('Flip',wnd,vbl+Param.RDK.FramesPerMove/NominalFrameRate-Slack);
        end
    end 
    
    
    %呈现刺激间隔  
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);    
    vbl = Screen('Flip',wnd, results(trial_i,15)+Param.RDK.Duration-Slack);
    results(trial_i,16)=vbl;

    %被试反应
    is_true = 0;
    while (is_true == 0 & GetSecs - results(trial_i,16) < Param.Trial.MaxRT)
        [keyIsDown_1,RT_time,keyCode] = KbCheck;
        if keyCode(Param.Keys.Right)
            results(trial_i,6) = 1;              %6 被试认为2为信号刺激
            if results(trial_i,5) > 0
                results(trial_i,8) = 1;          %7 被试反应正确与否
            end
            results(trial_i,9) = RT_time - results(trial_i,16);  % 9反应时
            is_true = 1;
        elseif keyCode(Param.Keys.Left)
            results(trial_i,6) = -1;              %被试认为1为信号刺激针
            if results(trial_i,4) > 0
                results(trial_i,8) = 1;
            end
            results(trial_i,9) = RT_time - results(trial_i,16);
            is_true = 1;
        elseif keyCode(Param.Keys.EscPress)
            abort;
        end
        %没有反应的试次会被删除
    end
    
    if ~results(trial_i,8)                  
        beep;
    end
    
    UD = PAL_AMUD_updateUD(UD,results(trial_i,8)); %update UD structure
    
    %呈现刺激间隔  
    Screen('Drawtext',wnd,'');   
    Screen('Flip',wnd);
    while (GetSecs - results(trial_i,10) < Param.Trial.Duration)
    end
    results(trial_i,17) = GetSecs-results(trial_i,10);
end
    
threshold_value = PAL_AMUD_analyzeUD(UD, 'reversals', max(UD.reversal)-Param.Staircase.ReversalsUsed);
disp('  ');
disp(['Threshold: ' num2str(threshold_value)]);
disp('  ');

%数据存储
CurrDir = pwd;
resultsDir = [CurrDir '\Results\DetectTrain\' SubjID '\'];
if ~isdir(resultsDir)                                                       %isder和mkder可以用来创建新文件夹
    mkdir(resultsDir);
end

cd(resultsDir);
results_name = [SubjID '_DetectTrain_results_session' num2str(session) '_block' num2str(block) '.mat'];
save(results_name,'results','UD','threshold_value');
cd(CurrDir); 

diff_temp = abs(results(:,4) - results(:,5));
trial_num_temp = max(results(:,1));
figure(1);
plot(1:trial_num_temp,diff_temp(1:trial_num_temp))

warning on;
%reset_test_gamma;
ShowCursor;
Screen('CloseAll');
 
delete *.asv  