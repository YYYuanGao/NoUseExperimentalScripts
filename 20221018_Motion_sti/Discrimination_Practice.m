function Discrimination_Practice(SubjID,session,block,MotionCoh)

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
%1��trial_i      %10���Դο�ʼʱ��=fixation����ʱ���
%2���̼�λ��       %11�������˶��̼���ʱ���
%3���˶�����       %12�������˶��̼���Ҫʱ��,ע�����ֵҪС��600ms
%4��ʵ�ʷ���1      %13���˶��̼�1����ʱ���
%5��ʵ�ʷ���2      %14���̼��������ʱ���
%6�����Է�Ӧ       %15���˶��̼�2����ʱ���
%7���ǶȲ���       %16��fixation����ʱ���
%8����ȷ���       %17���Դγ���
%9����Ӧʱ��       %18��coherence

results = zeros(Param.DisPrac.TrialNum,17);
trial_index = randperm(Param.DisPrac.TrialNum);
trial_index = mod(trial_index,Param.DisPrac.DirectionNum);

for trial_i = 1:Param.DisPrac.TrialNum
    
    %draw fixation
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);      
    vbl = Screen('Flip',wnd);
    results(trial_i,10) = vbl;
    
    results(trial_i,1) = trial_i;                             %1�Դ���
    Curr_cond = trial_index(trial_i)+1;   
    results(trial_i,2) = Param.DisPrac.StiLocation;           %2�̼�λ��
    results(trial_i,3) = Param.DisPrac.Directions(Curr_cond); %3�˶�����
    
    task_dif = Param.DisPrac.AngleDelta; 
    results(trial_i,7) = task_dif;
    results(trial_i,18) = MotionCoh;
    
    trial_order = ((rand-0.5)>=0);          % 0�ȳ��ֱ�׼�̼� 1����ֱ�׼�̼�
    trial_delta_ori= sign(rand-0.5);        % -1�Ǳ�׼�̼�����ڱ�׼�̼���ʱ����ת 1�Ǳ�׼�̼�����ڱ�׼�̼�˳ʱ����ת ��Ҫȷ�� 
    if trial_delta_ori == 0
        trial_delta_ori = 1;
    end

    if trial_order
        results(trial_i,5) = results(trial_i,3) + (rand-0.5)*2*Param.VarDir;  
        results(trial_i,4) = results(trial_i,5) + trial_delta_ori * task_dif;    
    else
        results(trial_i,4) = results(trial_i,3) + (rand-0.5)*2*Param.VarDir; 
        results(trial_i,5) = results(trial_i,4) + trial_delta_ori * task_dif;            
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
    
    if ~results(trial_i,8)                  
        beep;
    end
    
    %���ִ̼����  
    Screen('Drawtext',wnd,'');   
    Screen('Flip',wnd);
    while (GetSecs - results(trial_i,10) < Param.Trial.Duration)
    end
    results(trial_i,17) = GetSecs-results(trial_i,10);
end
    
prac_acc = sum(results(:,8))./Param.DisPrac.TrialNum;
disp('  ');
disp(['Accuracy: ' num2str(prac_acc)]);
disp('  ');

%���ݴ洢
CurrDir = pwd;
resultsDir = [CurrDir '\Results\DisPrac\' SubjID '\'];
if ~isdir(resultsDir)                                                       %isder��mkder���������������ļ���
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