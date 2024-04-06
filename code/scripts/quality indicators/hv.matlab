clear
clc

filepath1 = 'C:\Experiments on the Scalability\our algorithms\DTLZ7\1K\25\HV\HV\';
filepath2 = 'C:\Experiments on the Scalability\greedy\DTLZ7\1K\25\HV\HV\';

gsemo_arr = [];
greedy_arr = [];

for k = 1:30
	gsemo_mat = load([filepath1,'EXP2_1K_25_HV_',num2str(k),'.mat']);
	gsemo_num = gsemo_mat.RecordHV(54);
	gsemo_arr = [gsemo_arr;gsemo_num];
end

gsemo_std = std(gsemo_arr)
gsemo_mean = mean(gsemo_arr)

greedy_mat = load([filepath2,'Greedy_EXP2_1K_25_HV.mat']);
greedy_num = greedy_mat.res_hv(1)
greedy_arr = greedy_num*ones(30,1);

[p,h,stats] = ranksum(gsemo_arr, greedy_arr, 'alpha', 0.05, 'tail', 'both');

h
