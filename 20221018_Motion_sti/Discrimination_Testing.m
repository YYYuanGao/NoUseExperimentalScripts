function Discrimination_Testing(SubjID,session,block,MotionDirection,MotionCoh)

% subj_id = 'jiake'
% by Ke Jia
% Last modified 2021/10/27 18:18

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
%1��trial_i      %10�� trial start time = fixation presentation time point
%2���̼�λ��       %11�������˶��̼���ʱ���
%3���˶�����       %12�������˶��̼���Ҫʱ��,ע�����ֵҪС��600ms
%4��ʵ�ʷ���1      %13���˶��̼�1����ʱ���
%5��ʵ�ʷ���2      %14���̼��������ʱ���
%6�����Է�Ӧ       %15���˶��̼�2����ʱ���
%7���ǶȲ���       %16��fixation����ʱ���
%8����ȷ���       %17���Դγ���
%9����Ӧʱ��       %18��coherence

results = zeros(Param.Staircase.MaxTrial,18);

task_dif = Param.DisTest.AngleDelta;
UD = PAL_AMUD_setupUD('up',Param.Staircase.Up,'down',Param.Staircase.Down);
UD = PAL_AMUD_setupUD(UD,'StepSizeDown',Param.Staircase.StepSizeDown,'StepSizeUp',Param.Staircase.StepSizeUp,'stopcriterion',Param.Staircase.StopCriterion,'stoprule',Param.Staircase.StopRule,'startvalue',task_dif,...
        'xMax',Param.Staircase.xMax,'xMin',Param.Staircase.xMin,'truncate','yes');
    
while (~UD.stop & (size(UD.response)<Param.Staircase.MaxTrial))    
    trial_i = size(UD.response,2)+1;
    results(trial_i,1) = trial_i;                           %1�Դ���
    
    %draw fixation
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);      
    vbl = Screen('Flip',wnd);    
    results(trial_i,10) = vbl;
                                   
%    Curr_cond = trial_index(trial_i)+1;   
    results(trial_i,2) = Param.DisTest.StiLocation;           %2�̼�λ��
    results(trial_i,3) = MotionDirection;                     %3�˶�����    
    results(trial_i,7) = UD.xCurrent;
    results(trial_i,18) = MotionCoh;
    
    trial_order = ((rand-0.5)>=0);          % 0�ȳ��ֱ�׼�̼� 1����ֱ�׼�̼�
    trial_delta_ori= sign(rand-0.5);        % -1�Ǳ�׼�̼�����ڱ�׼�̼���ʱ����ת 1�Ǳ�׼�̼�����ڱ�׼�̼�˳ʱ����ת ��Ҫȷ�� 
    if trial_delta_ori == 0
        trial_delta_ori = 1;
    end

    if trial_order
        results(trial_i,5) = results(trial_i,3) + (rand-0.5)*2*Param.VarDir;  
        results(trial_i,4) = results(trial_i,5) + trial_delta_ori * UD.xCurrent;    
    else
        results(trial_i,4) = results(trial_i,3) + (rand-0.5)*2*Param.VarDir; 
        results(trial_i,5) = results(trial_i,4) + trial_delta_ori * UD.xCurrent;            
    end

    MyDot_1 = Gen_DotMatrix_bm(results(trial_i,4),results(trial_i,18),Param.RDK.NumFrames,Param.RDK.DotNum,Param.RDK.InnerRadius,Param.RDK.OuterRadius,Param.RDK.StepPerMove);
    MyDot_2 = Gen_DotMatrix_bm(results(trial_i,5),results(trial_i,18),Param.RDK.NumFrames,Param.RDK.DotNum,Param.RDK.InnerRadius,Param.RDK.OuterRadius,Param.RDK.StepPerMove);    
    results(trial_i,11) = GetSecs; 
    results(trial_i,12) = results(trial_i,11)-results(trial_i,10); 

    %�����˶���̼�
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
    
    
    %���ִ̼����  
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);    
    vbl = Screen('Flip',wnd, results(trial_i,13)+Param.RDK.Duration-Slack);
    results(trial_i,14)=vbl;
    
    
    %�����˶���̼�
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
    
    
    %���ִ̼����  
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);    
    vbl = Screen('Flip',wnd, results(trial_i,15)+Param.RDK.Duration-Slack);
    results(trial_i,16)=vbl;

    %���Է�Ӧ
    is_true = 0;
    while (is_true == 0 & GetSecs - results(trial_i,16) < Param.Trial.MaxRT)
        [keyIsDown_1,RT_time,keyCode] = KbCheck;
        if keyCode(Param.Keys.Right)
            results(trial_i,6) = 1;              %6 ������Ϊ2��1��˳ʱ��
            if results(trial_i,5)> results(trial_i,4)
                results(trial_i,8) = 1;          %7 ���Է�Ӧ��ȷ���
            end
            results(trial_i,9) = RT_time - results(trial_i,16);  % 9��Ӧʱ
            is_true = 1;
        elseif keyCode(Param.Keys.Left)
            results(trial_i,6) = -1;              %������Ϊ2��1����ʱ��
            if results(trial_i,5) < results(trial_i,4)
                results(trial_i,8) = 1;
            end
            results(trial_i,9) = RT_time - results(trial_i,16);
            is_true = 1;
        elseif keyCode(Param.Keys.EscPress)
            abort;
        end
        %û�з�Ӧ���Դλᱻɾ��
    end
    
    %if ~results(trial_i,8)                  
        %beep;
    %end
    
    UD = PAL_AMUD_updateUD(UD,results(trial_i,8)); %update UD structure
    
    %���ִ̼����  
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

%���ݴ洢
CurrDir = pwd;
resultsDir = [CurrDir '\Results\DisTest\' SubjID '\'];
if ~isdir(resultsDir)                                                       %isder��mkder���������������ļ���
    mkdir(resultsDir);
end

cd(resultsDir);
results_name = [SubjID '_DisTest_results_session' num2str(session) '_block' num2str(block) '.mat'];
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