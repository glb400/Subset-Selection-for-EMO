function [RecordHV,RecordR2,Subset,time]=GSEMOHV(PopObj,selNum,r,W)
    tic;
    profile on;
    %get size of Population.Objs
    % n = size(PopObj, 1);
    [n,m]=size(PopObj);
    %initialize the candidate solution set (called "population"): generate a Boolean string with all 0s (called "solution").
    population=zeros(1,n);
    %popSize: record the number of solutions in the population.
    popSize=1;
    %fitness: record the two objective values of a solution.
    fitness=zeros(1,2);
    %the first objective is f; for the special solution 00...00 (i.e., it does not select any variable) and the solutions with the number of selected variables not smaller than 2*k, set its first objective value as inf.  
    fitness(1)=0;
    %the second objective is the number of selected variables, i.e., the sum of all the bits.
    fitness(2)=sum(population); 

    % Record=[];
    RecordHV=[];
    RecordR2=[];

    % %generate direction matrix
    % num_vec=1000;
    % [W,~]=UniformVector(num_vec,m);




    r=r*max(PopObj,[],1);
    %repeat to improve the population; the number of iterations is set as 2*e*k^2*n suggested by our theoretical analysis.
    T=round(n*selNum*selNum*2*exp(1));
    % T=10*selNum*n;
    for t=1:T
        if mod(t,(selNum*n)) == 0
            disp('yes');
            temp=find(fitness(:,2)<=selNum);
            j=max(fitness(temp,2));
            seq=find(fitness(:,2)==j);
            selectedVariables=population(seq,:);
            %% output
            Subset=[];
            for i=1:n
                if selectedVariables(i) == 1
                    Subset=[Subset; PopObj(i,:)];
                end
            end
            % res=-1*fitness(seq,1);
            res1=HV(Subset,r);
            res2=R2ind(Subset,W,r);
            % res=Calcus(W,m,num_vec,Subset,j,r);
            disp(res1);
            disp(res2);
            % Record=[Record,res];
            RecordHV=[RecordHV,HV(Subset,r)];
            RecordR2=[RecordR2,R2ind(Subset,W,r)];
        end

        %randomly select a solution from the population and mutate it to generate a new solution.
        offspringChosen=population(randi([1,popSize],1,1),:);
        offspring=abs(offspringChosen-randsrc(1,n,[1,0; 1/n,1-1/n]));

        %compute the fitness of the new solution.
        offspringFit=zeros(1,2);
        offspringFit(2)=sum(offspring);

        if offspringFit(2)>=2*selNum
            continue;
        elseif offspringFit(2)==0
            offspringFit(1)=0;
        else
            offspringSet=find(offspring~=0);
            offspringNum=length(offspringSet);
            offspringElem=[];

            for i=1:offspringNum
                offspringElem=[offspringElem;PopObj(offspringSet(i),:)];
            end

            offspringFit(1)=R2ind(offspringElem,W,r);
        end

        %use the new solution to update the current population.
        if sum((fitness(1:popSize,1)>offspringFit(1)).*(fitness(1:popSize,2)<=offspringFit(2)))+sum((fitness(1:popSize,1)>=offspringFit(1)).*(fitness(1:popSize,2)<offspringFit(2)))>0
            continue;
        else
            deleteIndex=((fitness(1:popSize,1)<=offspringFit(1)).*(fitness(1:popSize,2)>=offspringFit(2)))'; 
        end
        %ndelete: record the index of the solutions to be kept.
        ndelete=find(deleteIndex==0);

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
    selectedVariables=population(seq,1:size(population,2));
    %% output
    Subset=[];
    for i=1:n
        if selectedVariables(i) == 1
            Subset=[Subset; PopObj(i,:)];
        end
    end
    time=toc;
    % res=-1*fitness(seq,1);
    res1=HV(Subset,r);
    res2=R2ind(Subset,W,r);
    % res=Calcus(W,m,num_vec,Subset,j,r);
    disp(res1);
    disp(res2);
end
