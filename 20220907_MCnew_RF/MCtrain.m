%% check parameters used
if exist([CurrDir '\Results\' SessName '\' SubjID '\' SubjID '_results_sess' num2str(Sess_Num) '_run' num2str(Run_Num) '.mat'],'file')
    disp(' ');
    disp([SubjID '_run' num2str(Run_Num) ' has been test, please enter a new run number!']); 
    disp(' ');
    abort;      
end

results = zeros(Param.RF.stoprule,15);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% results:
% 1 trial number              %10    trial onset
% 2 sti_loc                   %11    fix onset
% 3 sti_ori                   %12    ITI
% 4 SN of st1                 %13    st1
% 5 SN of st2                 %14    ISI
% 6 response                  %15    st2
% 7 task_diff
% 8 accuracy
% 9 reaction time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% quest settings
% if Sess_Num == 1 && Run_Num == 1
    task_dif = Param.RF.start;
% elseif Sess_Num == 1 && Run_Num ~= 1
%     previou_results = load([CurrDir '\Results\Train\' SubjID '\' SubjID '_train_results_sess' num2str(Sess_Num) '_run' num2str(Run_Num-1) '_ori' num2str(train_ori) '.mat']);   
% end

RF = PAL_AMRF_setupRF('priorAlphaRange', Param.RF.alphas, 'prior', Param.RF.prior,...
    'stopcriterion',Param.RF.stopcriterion,'stoprule',Param.RF.stoprule,'beta',Param.RF.beta,...
    'gamma',Param.RF.gamma,'lambda',Param.RF.lambda,'PF',Param.RF.PFfit,'meanmode',Param.RF.meanmode);
RF = PAL_AMRF_setupRF(RF,'startvalue',task_dif,'xMax',Param.RF.xMax,'xMin',Param.RF.xMin);

% start  
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

%% Go!
while (~RF.stop & (size(RF.response)<Param.RF.stoprule))
    trial_i = size(RF.response,2)+1;
    results(trial_i,1) = trial_i;
    results(trial_i,7) = RF.xCurrent;  
    results(trial_i,2) = location_used;
    results(trial_i,3) = Param.Stimuli.GratingOri(ori_used); 
   
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
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
    vbl = Screen('Flip',wnd);
    results(trial_i,11) = vbl-trial_onset; 
  
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
    Gabor_sti = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle)+x*cos(angle))+phase));
    Gabor_sti(sqrt(x.^2 + y.^2) > Param.Stimuli.OuterSize) = gray;
    
    Gabor_sti_final = repmat(Gabor_sti,[1,1,3]); 
    Gabor_sti_final(:,:,4) = annulus*white;
    
    mm = Screen('MakeTexture', wnd, Gabor_sti_final);     

    Screen('DrawTexture', wnd, mm,[],[Param.Stimuli.Locations(results(trial_i,2),1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)+Param.Stimuli.OuterSize]);
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
    vbl = Screen('Flip',wnd,vbl+Param.Trial.ITI); %-Slack
    results(trial_i,12) = vbl-trial_onset-results(trial_i,11); 

    %% ISI
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
    vbl = Screen('Flip',wnd,vbl+Param.Trial.StiDura); %-Slack
    results(trial_i,13) = vbl-trial_onset-sum(results(trial_i,11:12)); 

    %% sti2
    [x,y]=meshgrid(-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize,-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize);
    phase = rand*2*pi; 
    angle = results(trial_i,5)/180*pi;
    annulus=ones(Param.Stimuli.OuterSize*2+1);
    edge_control = sqrt(x.^2 + y.^2)./Param.Settings.PixelPerDegree;
    overrado=find(edge_control>Param.Stimuli.InnerRadius);

    len=length(overrado);
    for i=1:len
        annulus(overrado(i))=(annulus(overrado(i)).*exp(-1.*((((edge_control(overrado(i))-Param.Stimuli.InnerRadius)*Param.Settings.PixelPerDegree).^2)/(2.*(Param.Stimuli.SmoothSD.^2)))));
    end

    % sin grating
    Gabor_sti = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle)+x*cos(angle))+phase));
    Gabor_sti(sqrt(x.^2 + y.^2) > Param.Stimuli.OuterSize) = gray;
    
    Gabor_sti_final = repmat(Gabor_sti,[1,1,3]); 
    Gabor_sti_final(:,:,4) = annulus*white;
    
    mm = Screen('MakeTexture', wnd, Gabor_sti_final); 
    Screen('DrawTexture', wnd, mm,[],[Param.Stimuli.Locations(results(trial_i,2),1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)+Param.Stimuli.OuterSize]);
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
    vbl = Screen('Flip',wnd,vbl+Param.Trial.ISI); %-Slack
    stimulus_onset = vbl;
    results(trial_i,14) = vbl-trial_onset-sum(results(trial_i,11:13));

    %% response collection
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
    vbl = Screen('Flip',wnd,vbl+Param.Trial.StiDura);%-Slack
    results(trial_i,15) = vbl-trial_onset-sum(results(trial_i,11:14));      
    
    is_true = 0;
    while (is_true == 0 & GetSecs - stimulus_onset < Param.Trial.MaxRT)
        [ifkey,RT_time,keyCode] = KbCheck;
        if keyCode(Param.Keys.Right)
            results(trial_i,6) = 2;              %6 clockwise
            if results(trial_i,5)> results(trial_i,4)
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

    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    vbl = Screen('Flip',wnd);%-Slack
    
    if strcmp(SessName,'Prac') && (Run_Num<3) 
        if ~results(trial_i,8)                  
            beep;
        end
    end

    RF = PAL_AMRF_updateRF(RF,results(trial_i,7), results(trial_i,8));
    while (GetSecs - trial_onset < Param.Trial.Duration)
    end         
end
       
%% save the data
resultsDir = [CurrDir '\Results\' SessName '\' SubjID '\'];
if ~isdir(resultsDir)                                                       %isder和mkder可以用来创建新文件夹
    mkdir(resultsDir);
end

cd(resultsDir);
curr_threshold = RF.mean;
results_name = [SubjID '_results_sess' num2str(Sess_Num) '_run' num2str(Run_Num) '.mat'];
save(results_name,'results','RF','Param','curr_threshold');
cd(CurrDir); 

%% plot
figure(2);
trial_num_temp = 1:length(RF.x);
plot(trial_num_temp,RF.x,'k');
hold on;

plot(trial_num_temp(RF.response == 1),RF.x(RF.response == 1),'ko', 'MarkerFaceColor','k');
plot(trial_num_temp(RF.response == 0),RF.x(RF.response == 0),'ko', 'MarkerFaceColor','w');
set(gca,'FontSize',16);
% axis([0 max(trial_num_temp)+1 min(RF.x)-(max(RF.x)-min(RF.x))/10 max(RF.x)+(max(RF.x)-min(RF.x))/10]);
axis([0 max(trial_num_temp)+1 0 max(RF.x)+(max(RF.x)-min(RF.x))/10]);
% line([1 length(RF.x)], [Param.RF.trueParams(1) Param.RF.trueParams(1)],'linewidth', 2, 'linestyle', '--', 'color','k');
xlabel('Trials');
ylabel('Angle Difference');   
%%
reset_test_gamma;
ShowCursor;
Screen('CloseAll');
 
delete *.asv