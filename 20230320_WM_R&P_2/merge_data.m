
a=[];
List =dir('E:\WM_R&P_2\Results\train\oxy\oxy_results_sess1_run*.mat');
k =length(List);
for i=1:k
    file_name{i}=List(i).name;
    temp=importdata(file_name{i});
    a=[a;temp];
end
path_pri = ['marge.mat'];
save(path_pri,'a')