clear
clc

filepath1 = 'C:\Experiments on the Scalability\our algorithms\DTLZ7\1K\30\IGDp\';
filepath2 = 'C:\Experiments on the Scalability\greedy\DTLZ7\1K\30\IGDp\';

gsemo_arr = [];
greedy_arr = [];

for k = 1:30
	gsemo_mat = load([filepath1,'EXP2_1K_30_IGDp_',num2str(k),'.mat']);
	gsemo_num = gsemo_mat.Record(54);
	gsemo_arr = [gsemo_arr;gsemo_num];
end

gsemo_std = std(gsemo_arr)
gsemo_mean = mean(gsemo_arr)

greedy_mat = load([filepath2,'Greedy_EXP2_1K_30_IGDp.mat']);
greedy_num = greedy_mat.res(1)
greedy_arr = greedy_num*ones(30,1);

[p,h,stats] = ranksum(gsemo_arr, greedy_arr, 'alpha', 0.05, 'tail', 'both');

h
