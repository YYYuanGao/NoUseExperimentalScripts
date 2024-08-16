% Outer sti
% [x,y]=meshgrid(-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize,-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize);
% phase = rand*2*pi; 
% angle = 60/180*pi;
% sti2_final = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle)+x*cos(angle))+phase).*exp(-(x.^2 + y.^2)/2/(Param.Stimuli.SmoothSD.^2)));
% sti2_final(sqrt(x.^2 + y.^2) > Param.Stimuli.OuterSize) = gray; %circle mask
% 
% mm = Screen('MakeTexture', wnd, sti2_final);   
% location_temp = 1;
% Screen('DrawTexture', wnd, mm,[],[Param.Stimuli.Locations(location_temp,1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(location_temp,2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(location_temp,1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(location_temp,2)+Param.Stimuli.OuterSize]);
% 
[x,y]=meshgrid(-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize,-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize);
phase = rand*2*pi; 
angle = 105/180*pi; 
annulus=ones(Param.Stimuli.OuterSize*2+1);
edge_control = sqrt(x.^2 + y.^2)./Param.Settings.PixelPerDegree;
overrado=find(edge_control>Param.Stimuli.InnerRadius);

len=length(overrado);
for i=1:len
    annulus(overrado(i))=(annulus(overrado(i)).*exp(-1.*((((edge_control(overrado(i))-Param.Stimuli.InnerRadius)*Param.Settings.PixelPerDegree).^2)/(2.*(Param.Stimuli.SmoothSD.^2)))));
end

% sin grating
Gabor_sti = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle)+x*cos(angle))+phase).*exp(-(x.^2 + y.^2)/(2.*Param.Stimuli.SmoothSD.^2)));
Gabor_sti(sqrt(x.^2 + y.^2) > Param.Stimuli.OuterSize) = gray;

Gabor_sti_final = repmat(Gabor_sti,[1,1,3]); 
Gabor_sti_final(:,:,4) = annulus*white;

mm = Screen('MakeTexture', wnd, Gabor_sti_final); 
 
location_temp = 4;
Screen('DrawTexture', wnd, mm,[],[Param.Stimuli.Locations(location_temp,1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(location_temp,2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(location_temp,1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(location_temp,2)+Param.Stimuli.OuterSize]);


% fixation
Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
vbl = Screen('Flip',wnd);

     
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