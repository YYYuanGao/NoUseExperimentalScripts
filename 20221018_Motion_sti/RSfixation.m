function RSfixation

% by Ke Jia
% 2021/10/27 21:18
% 静息态实验使用的注视点刺激

parameter;


%*************************************
%*************Exp Start***************
%*************************************
Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
Screen('Flip', wnd);

is_true = 0;
while (is_true == 0)
    [~,~,keyCode] = KbCheck;
    if keyCode(Param.Keys.EscPress)
        abort;
    end
end

delete *.asv
 
        