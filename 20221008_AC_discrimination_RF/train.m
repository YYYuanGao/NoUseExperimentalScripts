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
% 1 trial number              %11    fixation
% 2 sti_loc                   %12    ITI duration
% 3 sti_ori                   %13    st1 duration
% 4 SN of st1                 %14    ISI
% 5 SN of st2                 %15    st2
% 6 response                  
% 7 task_diff                
% 8 accuracy                  
% 9 reaction time             
%10 trial onset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% quest settings
RF1 = PAL_AMRF_setupRF('priorAlphaRange', Param.RF.alphas,...
    'stopcriterion',Param.RF.stopcriterion,'stoprule',Param.RF.stoprule,'beta',Param.RF.beta,...
    'gamma',Param.RF.gamma,'lambda',Param.RF.lambda,'PF',Param.RF.PFfit,'meanmode',Param.RF.meanmode);
RF1 = PAL_AMRF_setupRF(RF1,'xMax',Param.RF.xMax,'xMin',Param.RF.xMin);

RF2 = PAL_AMRF_setupRF('priorAlphaRange', Param.RF.alphas,...
    'stopcriterion',Param.RF.stopcriterion,'stoprule',Param.RF.stoprule,'beta',Param.RF.beta,...
    'gamma',Param.RF.gamma,'lambda',Param.RF.lambda,'PF',Param.RF.PFfit,'meanmode',Param.RF.meanmode);
RF2 = PAL_AMRF_setupRF(RF2,'xMax',Param.RF.xMax,'xMin',Param.RF.xMin);

RF3 = PAL_AMRF_setupRF('priorAlphaRange', Param.RF.alphas,...
    'stopcriterion',Param.RF.stopcriterion,'stoprule',Param.RF.stoprule,'beta',Param.RF.beta,...
    'gamma',Param.RF.gamma,'lambda',Param.RF.lambda,'PF',Param.RF.PFfit,'meanmode',Param.RF.meanmode);
RF3 = PAL_AMRF_setupRF(RF3,'xMax',Param.RF.xMax,'xMin',Param.RF.xMin);

RF4 = PAL_AMRF_setupRF('priorAlphaRange', Param.RF.alphas,...
    'stopcriterion',Param.RF.stopcriterion,'stoprule',Param.RF.stoprule,'beta',Param.RF.beta,...
    'gamma',Param.RF.gamma,'lambda',Param.RF.lambda,'PF',Param.RF.PFfit,'meanmode',Param.RF.meanmode);
RF4 = PAL_AMRF_setupRF(RF4,'xMax',Param.RF.xMax,'xMin',Param.RF.xMin);

RF5 = PAL_AMRF_setupRF('priorAlphaRange', Param.RF.alphas,...
    'stopcriterion',Param.RF.stopcriterion,'stoprule',Param.RF.stoprule,'beta',Param.RF.beta,...
    'gamma',Param.RF.gamma,'lambda',Param.RF.lambda,'PF',Param.RF.PFfit,'meanmode',Param.RF.meanmode);
RF5 = PAL_AMRF_setupRF(RF5,'xMax',Param.RF.xMax,'xMin',Param.RF.xMin);

RF6 = PAL_AMRF_setupRF('priorAlphaRange', Param.RF.alphas,...
    'stopcriterion',Param.RF.stopcriterion,'stoprule',Param.RF.stoprule,'beta',Param.RF.beta,...
    'gamma',Param.RF.gamma,'lambda',Param.RF.lambda,'PF',Param.RF.PFfit,'meanmode',Param.RF.meanmode);
RF6 = PAL_AMRF_setupRF(RF6,'xMax',Param.RF.xMax,'xMin',Param.RF.xMin);


if Run_Num == 1 
    RF1 = PAL_AMRF_setupRF(RF1,'prior', Param.RF.prior,'stoprule',Param.RF.stoprule1);
    RF2 = PAL_AMRF_setupRF(RF2,'prior', Param.RF.prior,'stoprule',Param.RF.stoprule1);
    RF3 = PAL_AMRF_setupRF(RF3,'prior', Param.RF.prior,'stoprule',Param.RF.stoprule1);
    RF4 = PAL_AMRF_setupRF(RF4,'prior', Param.RF.prior,'stoprule',Param.RF.stoprule1);
    RF5 = PAL_AMRF_setupRF(RF5,'prior', Param.RF.prior,'stoprule',Param.RF.stoprule1);
    RF6 = PAL_AMRF_setupRF(RF6,'prior', Param.RF.prior,'stoprule',Param.RF.stoprule1);

elseif Run_Num ~= 1
    previou_results = load([CurrDir '\Results\' SessName '\' SubjID '\' SubjID '_results_sess' num2str(Sess_Num)  '_run' num2str(Run_Num-1) '.mat']);
    RF1_prior_temp = PAL_pdfNormal(Param.RF.alphas,previou_results.RF1.mean,Param.RF.sd);
    RF2_prior_temp = PAL_pdfNormal(Param.RF.alphas,previou_results.RF2.mean,Param.RF.sd);
    RF3_prior_temp = PAL_pdfNormal(Param.RF.alphas,previou_results.RF3.mean,Param.RF.sd);
    RF4_prior_temp = PAL_pdfNormal(Param.RF.alphas,previou_results.RF1.mean,Param.RF.sd);
    RF5_prior_temp = PAL_pdfNormal(Param.RF.alphas,previou_results.RF2.mean,Param.RF.sd);
    RF6_prior_temp = PAL_pdfNormal(Param.RF.alphas,previou_results.RF3.mean,Param.RF.sd);    
    RF1 = PAL_AMRF_setupRF(RF1,'prior', RF1_prior_temp,'stoprule',Param.RF.stoprule);
    RF2 = PAL_AMRF_setupRF(RF2,'prior', RF2_prior_temp,'stoprule',Param.RF.stoprule);
    RF3 = PAL_AMRF_setupRF(RF3,'prior', RF3_prior_temp,'stoprule',Param.RF.stoprule); 
    RF4 = PAL_AMRF_setupRF(RF4,'prior', RF4_prior_temp,'stoprule',Param.RF.stoprule);
    RF5 = PAL_AMRF_setupRF(RF5,'prior', RF5_prior_temp,'stoprule',Param.RF.stoprule);
    RF6 = PAL_AMRF_setupRF(RF6,'prior', RF6_prior_temp,'stoprule',Param.RF.stoprule);
end

%% start
DrawFormattedText(wnd,'Press space to start!','center','center', black);
Screen('Flip',wnd);
is_true = 0;
while (is_true == 0)
    [ifkey,RT_time,keyCode] = KbCheck;
    if keyCode(Param.Keys.Space)
        is_true = 1;
        sp.sendTrigger(Param.start);
    elseif keyCode(Param.Keys.Esc)
        abort;
    end
end
pause(0.5);

%% Go!
 while (~RF1.stop) || (~RF2.stop) || (~RF3.stop) || (~RF4.stop) || (~RF5.stop) || (~RF6.stop)
     trial_i = size(RF1.response,2)+size(RF2.response,2)+size(RF3.response,2)+size(RF4.response,2)+size(RF5.response,2)+size(RF6.response,2)+1;

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
     if (~RF4.stop)
         RF_continue = [RF_continue,4];
     end

     if (~RF5.stop)
         RF_continue = [RF_continue,5];
     end

     if (~RF6.stop)
         RF_continue = [RF_continue,6];
     end
     
     

    results(trial_i,1) = trial_i;
    results(trial_i,2) = location_used;
%     ori_temp = Shuffle(RF_continue);
%     ori_temp = ori_temp(1);
%     results(trial_i,3) = Param.Stimuli.GratingOri(ori_temp);

    if Sess_Num == 1
        ori_temp = Shuffle(RF_continue);
        ori_temp = ori_temp(1);
%     elseif Sess_Num == 2
%         ori_temp = mod((trial_i-1),3)+1;
    end
    results(trial_i,3) = Param.Stimuli.GratingOri(ori_temp);

    if ori_temp == 1
        results(trial_i,7) = RF1.xCurrent;
    elseif ori_temp == 2
        results(trial_i,7) = RF2.xCurrent;
    elseif  ori_temp == 3
        results(trial_i,7) = RF3.xCurrent;
    elseif  ori_temp == 4
        results(trial_i,7) = RF4.xCurrent;
    elseif  ori_temp == 5
        results(trial_i,7) = RF5.xCurrent;
    else
        results(trial_i,7) = RF6.xCurrent;
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

    %% ITI1
    rhyme = randperm(length(Param.Trial.ITI),1);
    results(trial_i,11) = Param.Trial.ITI(rhyme);
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
    vbl = Screen('Flip',wnd);
    sp.sendTrigger(Param.ITI1);
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
    Gabor_sti = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle)+x*cos(angle))+phase));
    Gabor_sti(sqrt(x.^2 + y.^2) > Param.Stimuli.OuterSize) = gray;
    
    Gabor_sti_final = repmat(Gabor_sti,[1,1,3]); 
    Gabor_sti_final(:,:,4) = annulus*white;
    
    mm = Screen('MakeTexture', wnd, Gabor_sti_final);     

    Screen('DrawTexture', wnd, mm,[],[Param.Stimuli.Locations(results(trial_i,2),1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)+Param.Stimuli.OuterSize]);
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
    vbl = Screen('Flip',wnd,vbl+Param.Trial.ITI); %-Slack
    if ori_temp == 1
    sp.sendTrigger(size(Param.Trial.sti1,1));
    elseif ori_temp == 2
    sp.sendTrigger(size(Param.Trial.sti1,2));
    elseif  ori_temp == 3
    sp.sendTrigger(size(Param.Trial.sti1,3));
    elseif  ori_temp == 4
    sp.sendTrigger(size(Param.Trial.sti1,4));
    elseif  ori_temp == 5
    sp.sendTrigger(size(Param.Trial.sti1,5));
    else
    sp.sendTrigger(size(Param.Trial.sti1,6));
    end
    results(trial_i,13) = vbl-trial_onset-sum(results(trial_i,12));  

    %% ISI
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
    vbl = Screen('Flip',wnd,vbl+Param.Trial.StiDura); %-Slack
    sp.sendTrigger(Param.ISI);
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
    sp.sendTrigger(Param.sti2);
    results(trial_i,15) = vbl-trial_onset-sum(results(trial_i,12:14));

    %% response collection
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
    vbl = Screen('Flip',wnd,vbl+Param.Trial.StiDura);%-Slack
    sp.sendTrigger(Param.response);
    results(trial_i,16) = vbl-trial_onset-sum(results(trial_i,12:15));      
    
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
                sp.sendTrigger(Param.correct);
            else     
                sp.sendTrigger(Param.wrong);
            end
            results(trial_i,9) = RT_time - stimulus_onset;
            is_true = 1;
        elseif keyCode(Param.Keys.Esc)
            abort;
        end
    end

    if ~results(trial_i,8)
                    beep;
%         DrawFormattedText(wnd,'Wrong!','center','center', black);        
%     else
%         DrawFormattedText(wnd,'Correct!','center','center', black);        
    end
    Screen('Flip',wnd);


    if ori_temp == 1
        RF1 = PAL_AMRF_updateRF(RF1,results(trial_i,7), results(trial_i,8));
    elseif ori_temp == 2
        RF2 = PAL_AMRF_updateRF(RF2,results(trial_i,7), results(trial_i,8));
    elseif ori_temp == 3
        RF3 = PAL_AMRF_updateRF(RF3,results(trial_i,7), results(trial_i,8));
    elseif ori_temp == 4
        RF4 = PAL_AMRF_updateRF(RF4,results(trial_i,7), results(trial_i,8));
    elseif ori_temp == 5
        RF5 = PAL_AMRF_updateRF(RF5,results(trial_i,7), results(trial_i,8));
    elseif ori_temp == 6
        RF6 = PAL_AMRF_updateRF(RF6,results(trial_i,7), results(trial_i,8));
    end
    while (GetSecs - trial_onset < Param.Trial.Duration)
    end      
   
    %% ITI2
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
    vbl = Screen('Flip',wnd);
    sp.sendTrigger(Param.ITI2);
    results(trial_i,17) = vbl-trial_onset-sum(results(trial_i,12:16));
end
       
%% save the data
resultsDir = [CurrDir '\Results\' SessName '\' SubjID '\'];
if ~isdir(resultsDir)                                                       %isder和mkder可以用来创建新文件夹
    mkdir(resultsDir);
end

cd(resultsDir);
curr_threshold1 = RF1.mean;
curr_threshold2 = RF2.mean;
curr_threshold3 = RF3.mean;
curr_threshold4 = RF4.mean;
curr_threshold5 = RF5.mean;
curr_threshold6 = RF6.mean;
results_name = [SubjID '_results_sess' num2str(Sess_Num) '_run' num2str(Run_Num) '.mat'];
save(results_name,'results','RF1','RF2','RF3','RF4','RF5','RF6','Param','curr_threshold1','curr_threshold2','curr_threshold3','curr_threshold4','curr_threshold5','curr_threshold6');
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
    elseif figure_i == 4
        trial_num_temp = 1:length(RF4.x);
        plot(trial_num_temp,RF4.x,'k');
        hold on;
        plot(trial_num_temp(RF4.response == 1),RF4.x(RF4.response == 1),'ko', 'MarkerFaceColor','k');
        plot(trial_num_temp(RF4.response == 0),RF4.x(RF4.response == 0),'ko', 'MarkerFaceColor','w');
        axis([0 max(trial_num_temp)+1 0 max(RF4.x)+(max(RF4.x)-min(RF4.x))/10]);
    elseif figure_i == 5
        trial_num_temp = 1:length(RF5.x);
        plot(trial_num_temp,RF5.x,'k');
        hold on;
        plot(trial_num_temp(RF5.response == 1),RF5.x(RF5.response == 1),'ko', 'MarkerFaceColor','k');
        plot(trial_num_temp(RF5.response == 0),RF5.x(RF5.response == 0),'ko', 'MarkerFaceColor','w');
        axis([0 max(trial_num_temp)+1 0 max(RF5.x)+(max(RF5.x)-min(RF5.x))/10]);
    elseif figure_i == 6
        trial_num_temp = 1:length(RF6.x);
        plot(trial_num_temp,RF6.x,'k');
        hold on;
        plot(trial_num_temp(RF6.response == 1),RF6.x(RF6.response == 1),'ko', 'MarkerFaceColor','k');
        plot(trial_num_temp(RF6.response == 0),RF6.x(RF6.response == 0),'ko', 'MarkerFaceColor','w');
        axis([0 max(trial_num_temp)+1 0 max(RF6.x)+(max(RF6.x)-min(RF6.x))/10]);

    end 
    set(gca,'FontSize',16);
    xlabel('Trials');
    ylabel('Angle Difference');  
end
%%
reset_test_gamma;
ShowCursor;
Screen('CloseAll');
 
delete *.asv