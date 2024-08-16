a = load('E:\20230320_WM_R&P_2\Results\train\czd\czd_results_sess1_run1.mat');
b = load('E:\20230320_WM_R&P_2\Results\train\jyl\jyl_results_sess1_run1.mat');
c = load('E:\20230320_WM_R&P_2\Results\train\qxy\qxy_results_sess1_run1.mat');
d = load('E:\20230320_WM_R&P_2\Results\train\zrj\zrj_results_sess1_run1.mat');
e = load('E:\20230320_WM_R&P_2\Results\train\zzf\zzf_results_sess1_run1.mat');



erge = [a.results; b.results;c.results;d.results;e.results];