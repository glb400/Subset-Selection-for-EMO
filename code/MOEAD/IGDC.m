function [IGDArray,Score] = IGDC(Population, PF, index, IGDArray)
    Subset = Population;
    Subset(index,:) = [];
    if size(Subset, 1) == 0
        Score = 1000 - IGD(Population, PF);
        IGDArray = min(pdist2(PF,Population),[],2);
    else
    Distance = pdist2(PF,Population(index, :));
    IGD_ori = mean(IGDArray);
    IGDArray = min(IGDArray, Distance);
    Score = IGD_ori - mean(IGDArray);
    end
end
    
   
               
