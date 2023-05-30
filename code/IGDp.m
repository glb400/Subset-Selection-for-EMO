function Score = IGDp(PopObj,PF)
    Distance = min(igddis(PF,PopObj),[],2);
    Score    = mean(Distance);
end

function [dis] = igddis(PF,Population)
    dis = zeros(size(PF, 1), size(Population, 1));
    for i = 1:size(PF, 1)
        for j = 1:size(Population, 1)
          dis(i, j) = sqrt(sum(max(Population(j,:) - PF(i,:), 0).^2));
        end
    end
end