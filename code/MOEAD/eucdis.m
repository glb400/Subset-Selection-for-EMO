function [dis] = eucdis(PF,Population)
    dis = zeros(size(PF, 1), size(Population, 1));
    for i = 1:size(PF, 1)
        for j = 1:size(Population, 1)
        dis(i, j) = sqrt(sum((Population(j,:) - PF(i,:)).^2));
        end
    end
end