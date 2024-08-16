%fixation_RS;

% wnd = Screen('OpenWindow',screenNumber,gray,Param.Settings.ScrnResolution); 
    Screen('FillRect',wnd,black,Param.Settings.ScrnResolution);
    Screen('Flip',wnd);

    pause(240);
    abort;
