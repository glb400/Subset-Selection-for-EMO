function [IGDArray,Score] = IGDCp(Population, PF, index, IGDArray)
    Subset = Population;
    Subset(index,:) = [];
    if size(Subset, 1) == 0
        Score = 1000 - IGDp(Population, PF);
        IGDArray = min(igddis(PF,Population),[],2);
    else
    Distance = igddis(PF,Population(index, :));
    IGD_ori = mean(IGDArray);
    IGDArray = min(IGDArray, Distance);
    Score = IGD_ori - mean(IGDArray);
    end
end

% function [dis] = igddis(PF,Population)
%     dis = zeros(size(PF, 1), size(Population, 1));
%     for i = 1:size(PF, 1)
%         for j = 1:size(Population, 1)
%         dis(i, j) = sqrt(sum(max(Population(j,:) - PF(i,:), 0).^2));
%         end
%     end
% end