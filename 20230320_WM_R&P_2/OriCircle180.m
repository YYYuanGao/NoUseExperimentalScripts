% inputori ~ [-180 360] 

% function outputori = OriCircle180(inputori)
% 
%     if inputori > 180
%         outputori = inputori - 180;
%     elseif inputori < 0
%         outputori = inputori + 180;
%     else
%         outputori = inputori;
%     end
%      
% end

function outputori = OriCircle180(inputori)
while inputori > 180 || inputori < 0
    if inputori > 180
        inputori = inputori - 180;
    elseif inputori < 0
        inputori = inputori + 180;
    end
end
outputori = inputori;
end