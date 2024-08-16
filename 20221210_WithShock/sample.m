clc;
clear;
HideCursor;
parameters;

% I1 = imread('Images2\A.png');
% I2=imread('o_1000006.bmp');
% b1 = Screen('MakeTexture',wnd,I1); %将图片I1读取到b1中备用
% b2=Screen('MakeTexture',winPt,I2);
% Screen('DrawTexture',wnd,b2); %将b1读取到刚刚新建的视窗中(还未显示)
% Screen('Flip',wnd); %切换为b1图片
% WaitSecs(2); 
% Screen('DrawTexture',wnd,b1);
% Screen('Flip',wnd);

load pic;
for i = 1:length(pic)
    pic{i,3} = imread(pic{i,2});
end

pic = repmat(pic,10,1);     %随机为原来的n倍
randIndex = randperm(length(pic));    %随机排序
pic = pic(randIndex,:);

%     j = 1:length(pic3)

j = randperm(length(pic),1)
file_image1 = Screen('MakeTexture',wnd,pic{j,3});  %将图片转为纹理数据
Screen('DrawTexture',wnd,file_image1);              %利用纹理数据绘制图片
Screen('Flip',wnd);

WaitSecs(1);

file_image1 = Screen('MakeTexture',wnd,pic{j,3});  %将图片转为纹理数据
Screen('DrawTexture',wnd,file_image1);              %利用纹理数据绘制图片
Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Jitter.Color, [], 1);
vbl = Screen('Flip',wnd);  %-Param.Screen.Slack


    %% close all
    is_true = 0;
    while ~is_true
        [~,~,keyCode] = KbCheck;
        if keyCode(Param.Keys.Esc)
            is_true = 1;
        end
    end


Screen('CloseAll');
reset_test_gamma;