classdef DRUPAL < PROBLEM
    properties(Access = public)
        drupals;  % drupal dataset
        Whv;  % direction vectors of hv approximation 
        Wr2;  % weight vectors of r2
        k;  % the size of subset
        r;  % the reference coefficient
    end
    methods
    	%% Default settings of the problem
        function Setting(obj)
            % Load data
            obj.drupals = importdata('C:\Users\yons\Desktop\NSGAII\Dataset\normalizedNondominatedSolution.mat');
            obj.Whv = importdata('C:\Users\yons\Desktop\NSGAII\REALWORLD\W_HV_REALWORLD.mat');
            obj.Wr2 = importdata('C:\Users\yons\Desktop\NSGAII\REALWORLD\W_R2_REALWORLD.mat');

            % Parameter setting
            obj.M = 2;
            obj.D = size(obj.drupals,1);
            obj.k = 10;
            obj.r = 1.1;

            obj.lower    = zeros(1, obj.D);
            obj.upper    = ones(1, obj.D);
            
            obj.maxFE    = 233300;
            obj.encoding = 'binary';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
			PopObj = zeros(size(PopDec,1),obj.M);

			% objective 1 is the size
			PopObj(:,1) = sum(PopDec,2);

			% objective 2 is the quality indicator value
			% % HV
			% for i = 1:length(PopObj(:,2))
			% 	elem = PopDec(i,:);
			% 	% get subset from drupal dataset
			% 	Subset = [];
			% 	for j = 1:obj.D
			% 		if elem(j) == 1
			% 			Subset = [Subset; obj.drupals(j,:)];
			% 		end
			% 	end
			% 	% set size = 0 or size >= 2k to be a number larger than supremum
			% 	if length(Subset) == 0 || length(Subset) >= 2 * obj.k
			% 		PopObj(i,2) = 0;
			% 	else
			% 		PopObj(i,2) = -1 * R2ind(Subset,obj.Whv,obj.r);
			% 		% PopObj(i,2) = -1 * HV(Subset,obj.r);	
			% 	end
			% end

			% % IGD
			% for i = 1:length(PopObj(:,2))
			% 	elem = PopDec(i,:);
			% 	% get subset from drupal dataset
			% 	Subset = [];
			% 	for j = 1:obj.D
			% 		if elem(j) == 1
			% 			Subset = [Subset; obj.drupals(j,:)];
			% 		end
			% 	end
			% 	% set size = 0 or size >= 2k to be a number larger than supremum
			% 	if length(Subset) == 0 || length(Subset) >= 2 * obj.k
			% 		PopObj(i,2) = 100;
			% 	else
			% 		PopObj(i,2) = IGD(Subset,obj.drupals);
			% 	end
			% end

			% % IGDp
			% for i = 1:length(PopObj(:,2))
			% 	elem = PopDec(i,:);
			% 	% get subset from drupal dataset
			% 	Subset = [];
			% 	for j = 1:obj.D
			% 		if elem(j) == 1
			% 			Subset = [Subset; obj.drupals(j,:)];
			% 		end
			% 	end
			% 	% set size = 0 or size >= 2k to be a number larger than supremum
			% 	if length(Subset) == 0 || length(Subset) >= 2 * obj.k
			% 		PopObj(i,2) = 100;
			% 	else
			% 		PopObj(i,2) = IGDp(Subset,obj.drupals);
			% 	end
			% end

			% R2
			for i = 1:length(PopObj(:,2))
				elem = PopDec(i,:);
				% get subset from drupal dataset
				Subset = [];
				for j = 1:obj.D
					if elem(j) == 1
						Subset = [Subset; obj.drupals(j,:)];
					end
				end
				% set size = 0 or size >= 2k to be a number larger than supremum
				if length(Subset) == 0 || length(Subset) >= 2 * obj.k
					PopObj(i,2) = 100;
				else
					PopObj(i,2) = R2Tchebycheff(Subset,obj.Wr2,obj.r);
				end
			end	
        end
    end
end
