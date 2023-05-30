function [nondominatedSolution]= deleteSolutions(data_set)
    nondominatedSolution = [];
    for i=1:length(data_set)
        flag = 1;
        di = data_set(i,:);
        for j=1:length(data_set)
            if i == j
                continue;
            end
            dj = data_set(j,:);
            if sum(di == dj) == 9
                continue;
            end
            compare1 = di(1:7) <= dj(1:7);
            compare2 = di(8:9) >= dj(8:9);
            if sum([compare1, compare2]) == 9
                flag = 0;
                break;
            end
        end
        if flag == 1
            nondominatedSolution = [nondominatedSolution; di];
        end
    end
    nondominatedSolution = nondominatedSolution(all(~isnan(nondominatedSolution),2),:);
end