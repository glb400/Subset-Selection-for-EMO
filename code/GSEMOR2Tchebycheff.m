function [RecordR2Tch,Subset,time]=GSEMOR2Tchebycheff(PopObj,selNum,r,W)
    tic;
    profile on;
    %get size of Population.Objs
    [n,m]=size(PopObj);
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


    RecordR2Tch=[];

    num_vec=size(W,1);
    % %generate weight matrix
    % [W,~]=UniformVector2(num_vec,m);
    
    % utopian point
    r=r*max(PopObj,[],1);





    %add 0-solution as initial population
    
    % IGDArray=(1./zeros(1,n))';
    gteArray=(1./zeros(1,num_vec))';



    %id set is too big if we use dictionary to store array Dw 
    %instead we add id into solution
    population = [population, 1];
    
    gteArraySet=(gteArray);
    % IGDArraySet = (IGDArray);


    % Record=[];





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


            [res,~]=R2Tchebycheff(Subset,W,r);
            disp(res);
            RecordR2Tch=[RecordR2Tch,res];


            

            % res=IGD(Subset,PopObj);
            % disp(res);
            % Record=[Record,res];




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



            % %get array D
            % IGDArrayId=PopObjChosen(:,size(population,2));
            % IGDArray=IGDArraySet(:,IGDArrayId);



            %get array Dw
            gteArrayId=PopObjChosen(:,size(population,2));
            gteArray=gteArraySet(:,gteArrayId);


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
                % IGDArray=(1./zeros(1,n))';

                gteArray=(1./zeros(1,num_vec))';




            elseif offspringChosenNum==0
                % IGDArray=min(eucdis(PopObj,offspringSet),[],2);

                [~,gteArray]=R2Tchebycheff(offspringSet,W,r);




            else
                %add calculation
                addSet=find((offspring-offspringChosen)==1);
                for i=1:length(addSet)
                    si=PopObj(addSet(i),:);
                    % offspringChosenSet=[offspringChosenSet;si];
                    [~,gteArraySi]=R2Tchebycheff(si,W,r);
                    gteArray=min([gteArray,gteArraySi],[],2);
                    % [IGDArray,~]=IGDC(offspringChosenSet, PopObj, size(offspringChosenSet,1), IGDArray);
                end

                %delete calculation
                deleteSet=find((offspring-offspringChosen)==-1);
                for i=1:length(deleteSet)
                    % disp('hi');

                    % %calculate Array D' with respect to deleted element s_{i}
                    % si=PopObj(deleteSet(i),:);


                    % Distance=[];



                    % for j=1:n
                    %     Distance=[Distance;eucdis(PopObj(j,:),si)];
                    % end
                    % recalculateSet=find((Distance-IGDArray)==0);




                    % for j=1:length(recalculateSet)
                    %     id=recalculateSet(j);
                    %     IGDArray(id)=inf;
                    %     for k=1:offspringNum
                    %         IGDArray(id)=min(IGDArray(id),eucdis(PopObj(id,:),PopObj(offspringElem(k),:)));
                    %     end
                    % end

                    %calculate Array D' with respect to deleted element s_{i}
                    si=PopObj(deleteSet(i),:);
                    [~, gteArraySi]=R2Tchebycheff(si,W,r);

                    % disp(size(gteArraySi));
                    % disp('---');
                    % disp(size(gteArray));

                    recalculateSet=find((gteArraySi-gteArray)==0);
                    
                    for j=1:length(recalculateSet)
                        % disp('hello');
                        id=recalculateSet(j);
                        gteArray(id)=inf;
                        temp = abs(offspringSet-r).*W(id,:);
                        gteArray(id) = min(max(temp,[],2));

                        % IGDArray(id)=inf;
                        % for k=1:offspringNum
                        %     IGDArray(id)=min(IGDArray(id),eucdis(PopObj(id,:),PopObj(offspringElem(k),:)));
                        % end
                    end




                end



            end
            % offspringFit(1)=mean(IGDArray);
            offspringFit(1)=mean(gteArray);
            % disp(offspringFit(1));

        end




        %use the new solution to update the current population.
        if sum((fitness(1:popSize,1)<offspringFit(1)).*(fitness(1:popSize,2)<=offspringFit(2)))+sum((fitness(1:popSize,1)<=offspringFit(1)).*(fitness(1:popSize,2)<offspringFit(2)))>0
            continue;
        else
            deleteIndex=((fitness(1:popSize,1)>=offspringFit(1)).*(fitness(1:popSize,2)>=offspringFit(2)))'; 
        end
        %ndelete: record the index of the solutions to be kept.
        ndelete=find(deleteIndex==0);




        % newIGDArraySet=[];
        % for i=1:length(ndelete)
        %     curId=population(ndelete(i),size(population,2));
        %     newIGDArraySet=[newIGDArraySet,IGDArraySet(:,curId)];
        %     population(ndelete(i),size(population,2))=i;
        % end
        % IGDArraySet=[newIGDArraySet,IGDArray];


        newgteArraySet=[];
        for i=1:length(ndelete)
            curId=population(ndelete(i),size(population,2));            
            newgteArraySet=[newgteArraySet,gteArraySet(:,curId)];
            population(ndelete(i),size(population,2))=i;
        end
        gteArraySet=[newgteArraySet,gteArray];


        % offspring=[offspring,size(IGDArraySet,2)]; 
        % population=[population(ndelete,:)',offspring']';
        % fitness=[fitness(ndelete,:)',offspringFit']';
        % popSize=size(population,1);

        offspring=[offspring,size(gteArraySet,2)]; 
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



    % res=IGD(Subset,PopObj);
    % disp(res);


    [res,~]=R2Tchebycheff(Subset,W,r);
    disp(res);

end
