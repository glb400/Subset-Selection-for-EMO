classdef NSGAII < ALGORITHM
% <multi> <real/binary/permutation> <constrained/none>
% Nondominated sorting genetic algorithm II

%------------------------------- Reference --------------------------------
% K. Deb, A. Pratap, S. Agarwal, and T. Meyarivan, A fast and elitist
% multiobjective genetic algorithm: NSGA-II, IEEE Transactions on
% Evolutionary Computation, 2002, 6(2): 182-197.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();
            [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
                Offspring  = OperatorGA(Population(MatingPool));
                [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Problem.N);


                % % output the optimum
                % % HV
                % bestValue = 0;
                % hvValue = 0;
                % objValues = Population.objs;
                % objValues(:,2) = -1 * objValues(:,2);
                % decValues = Population.decs;
                % for i = 1:Problem.N
                %     if objValues(i,1) <= 10
                %         if bestValue < objValues(i,2)
                %             bestValue = objValues(i,2);
                %             % get subset from drupal dataset
                %             elem = decValues(i,:);
                %             Subset = [];
                %             for j = 1:Problem.D
                %                 if elem(j) == 1
                %                     Subset = [Subset; Problem.drupals(j,:)];
                %                 end
                %             end
                %             hvValue = HV(Subset,1.1*ones(1,9));
                %         end
                %     end
                % end
                % disp(bestValue);
                % save('NSGAII_EXP3_HV_30.mat', 'hvValue');
                % save('NSGAII_EXP3_R2APPROX_30.mat', 'bestValue');

                % % IGD
                % bestValue = inf;
                % objValues = Population.objs;
                % for i = 1:Problem.N
                %     if objValues(i,1) <= 10
                %         if bestValue > objValues(i,2)
                %             bestValue = objValues(i,2);
                %         end
                %     end
                % end
                % disp(bestValue);
                % save('NSGAII_EXP3_IGD_30.mat', 'bestValue');

                % % IGDp
                % bestValue = inf;
                % objValues = Population.objs;
                % for i = 1:Problem.N
                %     if objValues(i,1) <= 10
                %         if bestValue > objValues(i,2)
                %             bestValue = objValues(i,2);
                %         end
                %     end
                % end
                % disp(bestValue);
                % save('NSGAII_EXP3_IGDp_30.mat', 'bestValue');

                % R2
                bestValue = inf;
                objValues = Population.objs;
                for i = 1:Problem.N
                    if objValues(i,1) <= 10
                        if bestValue > objValues(i,2)
                            bestValue = objValues(i,2);
                        end
                    end
                end
                disp(bestValue);
                save('NSGAII_EXP3_R2_30.mat', 'bestValue');
            end
        end
    end
end