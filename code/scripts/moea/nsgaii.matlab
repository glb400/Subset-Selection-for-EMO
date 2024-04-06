clear
clc

for k = 1:30
	rng('shuffle')
	A = DRUPAL();
	% run moead
	platemo('algorithm',@NSGAII,'problem',@DRUPAL,'save',10);
end

