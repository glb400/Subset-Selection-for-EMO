function [Record,Subset,time]=GSEMOIGD(PopObj,selNum)
%%  ATTN
%   ATTN: This package is free for academic usage. You can run it at your own risk. For other purposes, please contact Prof. Zhi-Hua Zhou (zhouzh@nju.edu.cn)
%%  ATTN2
%   ATTN2: This package was developed by Mr. Chao Qian (qianc@lamda.nju.edu.cn). For any problem concerning the code, please feel free to contact Mr. Qian.
%%  Some varables used in the code
%   input: 
%      NOTE THAT we use a Boolean string "1001....1100" to represent a subset of variables, where "1" means that the corresponding variable is selected, and "0" means ignoring the variable.
%      n: the total number of variables. 
%      f: a given criterion to be optimized; its input should be a subset of variables, i.e., a Boolean string of length n; its output should be a real value.
%      k: the constraint on the number of selected variables.
%      NOTE THAT we assume that f is to be minimized.
%   ouptut: 
%      selectedVariables: a Boolean string of length n representing the selected variables, the number of which is not larger than k. 

    tic;
    profile on;
    %get size of Population.Objs
    n = size(PopObj, 1);
    %initialize the candidate solution set (called "population"): generate a Boolean string with all 0s (called "solution").
    population=zeros(1,n);
    %popSize: record the number of solutions in the population.
    popSize=1;
    %fitness: record the two objective values of a solution.
    fitness=zeros(1,2);
    %the first objective is f; for the special solution 00...00 (i.e., it does not select any variable) and the solutions with the number of selected variables not smaller than 2*k, set its first objective value as inf.  
    fitness(1)=inf;
    %the second objective is the number of selected variables, i.e., the sum of all the bits.
    fitness(2)=sum(population); 

    %add 0-solution as initial population
    IGDArray=(1./zeros(1,n))';

    %id set is too big if we use dictionary to store array D
    % IGDArrayDict = struct(num2str(population),IGDArray);
    %instead we add id into solution
    population = [population, 1];
    IGDArraySet = (IGDArray);

    Record=[];

    %repeat to improve the population; the number of iterations is set as 2*e*k^2*n suggested by our theoretical analysis.
    T=round(n*selNum*selNum*2*exp(1));
    % T=10*selNum*n;
    for t=1:T
        if mod(t,(selNum*n)) == 0
            disp('yes');
            temp=find(fitness(:,2)<=selNum);
            j=max(fitness(temp,2));
            seq=find(fitness(:,2)==j);
            selectedVariables=population(seq,1:size(population,2)-1);
            %% output
            Subset=[];
            for i=1:n
                if selectedVariables(i) == 1
                    Subset=[Subset; PopObj(i,:)];
                end
            end
            res=IGD(Subset,PopObj);
            disp(res);
            Record=[Record,res];
        end


        %randomly select a solution from the population and mutate it to generate a new solution.
        PopObjChosen=population(randi([1,popSize],1,1),:);

        offspringChosen=PopObjChosen(1:size(population,2)-1);
        offspring=abs(offspringChosen-randsrc(1,n,[1,0; 1/n,1-1/n]));

        %compute the fitness of the new solution.
        offspringFit=zeros(1,2);
        offspringFit(2)=sum(offspring);
        
        if offspringFit(2)>=2*selNum
            continue;
        elseif offspringFit(2)==0
            offspringFit(1)=inf;
        else
            %get array D
            IGDArrayId=PopObjChosen(:,size(population,2));
            IGDArray=IGDArraySet(:,IGDArrayId);

            %get offspring/offspringChosen information
            offspringElem=find(offspring~=0);
            offspringChosenElem=find(offspringChosen~=0);
            offspringNum=length(offspringElem);
            offspringChosenNum=length(offspringChosenElem);

            %get offspringSet/offspringChosenSet/offspringTotalSet
            offspringSet=[];
            offspringChosenSet=[];
            for i=1:offspringNum
                offspringSet=[offspringSet;PopObj(offspringElem(i),:)];
            end
            for i=1:offspringChosenNum
                offspringChosenSet=[offspringChosenSet;PopObj(offspringChosenElem(i),:)];
            end       

            if offspringNum==0
                IGDArray=(1./zeros(1,n))';
            elseif offspringChosenNum==0
                IGDArray=min(eucdis(PopObj,offspringSet),[],2);
            else
                %add calculation
                addSet=find((offspring-offspringChosen)==1);
                for i=1:length(addSet)
                    si=PopObj(addSet(i),:);
                    offspringChosenSet=[offspringChosenSet;si];
                    [IGDArray,~]=IGDC(offspringChosenSet, PopObj, size(offspringChosenSet,1), IGDArray);
                end

                %delete calculation
                deleteSet=find((offspring-offspringChosen)==-1);
                for i=1:length(deleteSet)
                    %calculate Array D' with respect to deleted element s_{i}
                    si=PopObj(deleteSet(i),:);
                    Distance=[];
                    for j=1:n
                        Distance=[Distance;eucdis(PopObj(j,:),si)];
                    end
                    recalculateSet=find((Distance-IGDArray)==0);

                    for j=1:length(recalculateSet)
                        id=recalculateSet(j);
                        IGDArray(id)=inf;
                        for k=1:offspringNum
                            IGDArray(id)=min(IGDArray(id),eucdis(PopObj(id,:),PopObj(offspringElem(k),:)));
                        end
                    end
                end
            end
            offspringFit(1)=mean(IGDArray);
        end

        %use the new solution to update the current population.
        if sum((fitness(1:popSize,1)<offspringFit(1)).*(fitness(1:popSize,2)<=offspringFit(2)))+sum((fitness(1:popSize,1)<=offspringFit(1)).*(fitness(1:popSize,2)<offspringFit(2)))>0
            continue;
        else
            deleteIndex=((fitness(1:popSize,1)>=offspringFit(1)).*(fitness(1:popSize,2)>=offspringFit(2)))'; 
        end
        %ndelete: record the index of the solutions to be kept.
        ndelete=find(deleteIndex==0);
        newIGDArraySet=[];
        for i=1:length(ndelete)
            curId=population(ndelete(i),size(population,2));
            newIGDArraySet=[newIGDArraySet,IGDArraySet(:,curId)];
            population(ndelete(i),size(population,2))=i;
        end
        IGDArraySet=[newIGDArraySet,IGDArray];


        offspring=[offspring,size(IGDArraySet,2)]; 
        population=[population(ndelete,:)',offspring']';
        fitness=[fitness(ndelete,:)',offspringFit']';
        popSize=size(population,1);
    end
    
    %select the final solution according to the constraint k on the number of selected variables.
    for i=1:popSize
        disp(fitness(i,1));
        disp(fitness(i,2));
        disp('---');
    end

    temp=find(fitness(:,2)<=selNum);
    j=max(fitness(temp,2));
    seq=find(fitness(:,2)==j);
    selectedVariables=population(seq,1:size(population,2)-1);
    
    %% output
    Subset=[];
    for i=1:n
        if selectedVariables(i) == 1
            Subset=[Subset; PopObj(i,:)];
        end
    end
    time=toc;
    res=IGD(Subset,PopObj);
    disp(res);
end
