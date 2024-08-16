load pic;
for i = 1:length(pic)
    pic{i,3} = imread(pic{i,2});
    pic{i,3}(pic{i,3}==0) = 128;
end

