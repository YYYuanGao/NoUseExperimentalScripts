pic = repmat(pic,100,1);     %随机为原来的n倍
randIndex = randperm(length(pic));    %随机排序
pic = pic(randIndex,:);