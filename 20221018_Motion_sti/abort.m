%终止实验程序
function abort  
    ShowCursor;              %显示鼠标
    Screen('CloseAll');      %关闭所有屏
    %reset_test_gamma;       %去除gamma矫正
    warning on;              %在命令窗口显示warning
    error('aborting due to user cancellation...');
end