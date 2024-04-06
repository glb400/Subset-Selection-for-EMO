clear
clc

filepath1 = 'C:\Experiments on Benchmark and Real-world Problems\our algorithms\EMOA\R2\M3\';
filepath2 = 'C:\Experiments on Benchmark and Real-world Problems\greedy\EMOA\R2\M3\';

greedy_arr = [];
gsemo_value_arr = [];

for i = 1:54	
	gsemo_arr = [];
	for k = 1:30
		gsemo_mat = load([filepath1,'BENCHMARK_NSGAII_R2_M3_minus_dtlz2_',num2str(k),'.mat']);
		gsemo_num = gsemo_mat.RecordR2Tch(i);
		gsemo_arr = [gsemo_arr;gsemo_num];
	end
	gsemo_value_arr = [gsemo_value_arr; mean(gsemo_arr)];
end

% gsemo_std = std(gsemo_arr)
% gsemo_mean = mean(gsemo_arr)

greedy_mat = load([filepath2,'Greedy_BENCHMARK_NSGAII_R2_M3_minus_dtlz2.mat']);
greedy_num = greedy_mat.res(1)
greedy_arr = greedy_num*ones(54,1);

figure;

plot(gsemo_value_arr, 'LineWidth', 2, 'Color', [1,0,0])
xlim = get(gca, 'Xlim');

hold on
plot(xlim, [greedy_num, greedy_num], 'k-', 'LineWidth', 2, 'Color', [0,0,1], 'linestyle', '--')
axis([1 54 0.650 0.660])
ylabel('R2 value')
xlabel('Running time in \sl{kn}')
legend('GSEMO-ACC', 'Greedy', 'Location', 'SouthEast')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
