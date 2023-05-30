classdef MOEAD < ALGORITHM
% <multi/many> <real/binary/permutation>
% Multiobjective evolutionary algorithm based on decomposition
% type --- 1 --- The type of aggregation function

%------------------------------- Reference --------------------------------
% Q. Zhang and H. Li, MOEA/D: A multiobjective evolutionary algorithm based
% on decomposition, IEEE Transactions on Evolutionary Computation, 2007,
% 11(6): 712-731.
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
            %% Parameter setting
            type = Algorithm.ParameterSet(1);

            %% Generate the weight vectors
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            T = ceil(Problem.N/10);

            %% Detect the neighbours of each solution
            B = pdist2(W,W);
            [~,B] = sort(B,2);
            B = B(:,1:T);

            %% Generate random population
            Population = Problem.Initialization();
            Z = min(Population.objs,[],1);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                % For each solution
                for i = 1 : Problem.N
                    % Choose the parents
                    P = B(i,randperm(size(B,2)));

                    % Generate an offspring
                    Offspring = OperatorGAhalf(Population(P(1:2)));

                    % Update the ideal point
                    Z = min(Z,Offspring.obj);

                    % Update the neighbours
                    switch type
                        case 1
                            % PBI approach
                            normW   = sqrt(sum(W(P,:).^2,2));
                            normP   = sqrt(sum((Population(P).objs-repmat(Z,T,1)).^2,2));
                            normO   = sqrt(sum((Offspring.obj-Z).^2,2));
                            CosineP = sum((Population(P).objs-repmat(Z,T,1)).*W(P,:),2)./normW./normP;
                            CosineO = sum(repmat(Offspring.obj-Z,T,1).*W(P,:),2)./normW./normO;
                            g_old   = normP.*CosineP + 5*normP.*sqrt(1-CosineP.^2);
                            g_new   = normO.*CosineO + 5*normO.*sqrt(1-CosineO.^2);
                        case 2
                            % Tchebycheff approach
                            g_old = max(abs(Population(P).objs-repmat(Z,T,1)).*W(P,:),[],2);
                            g_new = max(repmat(abs(Offspring.obj-Z),T,1).*W(P,:),[],2);
                        case 3
                            % Tchebycheff approach with normalization
                            Zmax  = max(Population.objs,[],1);
                            g_old = max(abs(Population(P).objs-repmat(Z,T,1))./repmat(Zmax-Z,T,1).*W(P,:),[],2);
                            g_new = max(repmat(abs(Offspring.obj-Z)./(Zmax-Z),T,1).*W(P,:),[],2);
                        case 4
                            % Modified Tchebycheff approach
                            g_old = max(abs(Population(P).objs-repmat(Z,T,1))./W(P,:),[],2);
                            g_new = max(repmat(abs(Offspring.obj-Z),T,1)./W(P,:),[],2);
                    end
                    Population(P(g_old>=g_new)) = Offspring;
                end

                % output the optimum
                % HV
                bestValue = inf;
                hvValue = 0;
                objValues = Population.objs;
                
                decValues = Population.decs;
                for i = 1:Problem.N
                    if objValues(i,1) <= (10 / 429) + 1e-6
                        if bestValue > objValues(i,2)
                            bestValue = objValues(i,2);
                            % get subset from drupal dataset
                            elem = decValues(i,:);
                            Subset = [];
                            for j = 1:Problem.D
                                if elem(j) == 1
                                    Subset = [Subset; Problem.drupals(j,:)];
                                end
                            end
                            hvValue = HV(Subset,1.1*ones(1,9));
                        end
                    end
                end
                disp(bestValue);
                save('MOEAD_EXP3_HV_30.mat', 'hvValue');
                save('MOEAD_EXP3_R2APPROX_30.mat', 'bestValue');

                % % IGD
                % bestValue = inf;
                % objValues = Population.objs;
                % for i = 1:Problem.N
                %     if objValues(i,1) <= (10 / 429) + 1e-6
                %         if bestValue > objValues(i,2)
                %             bestValue = objValues(i,2);
                %         end
                %     end
                % end
                % disp(bestValue);
                % save('MOEAD_EXP3_IGD_30.mat', 'bestValue');

                % % IGDp
                % bestValue = inf;
                % objValues = Population.objs;
                % for i = 1:Problem.N
                %     if objValues(i,1) <= (10 / 429) + 1e-6
                %         if bestValue > objValues(i,2)
                %             bestValue = objValues(i,2);
                %         end
                %     end
                % end
                % disp(bestValue);
                % save('MOEAD_EXP3_IGDp_30.mat', 'bestValue');

                % % R2
                % bestValue = inf;
                % objValues = Population.objs;
                % for i = 1:Problem.N
                %     if objValues(i,1) <= (10 / 429) + 1e-6
                %         if bestValue > objValues(i,2)
                %             bestValue = objValues(i,2);
                %         end
                %     end
                % end
                % disp(bestValue);
                % save('MOEAD_EXP3_R2_30.mat', 'bestValue');                      
            end
        end
    end
end