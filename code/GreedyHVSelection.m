function [Subset,res_hv,res_r2] = GreedyHVSelection(PopObj,k,r,W)
    a = size(PopObj, 1);

    selectedPop = [];
    currentHV = -inf;
    ids = zeros(1,a);

    r = r * max(PopObj, [], 1);

    for i = 1 : k
        elemId = 0;
        for n = 1 : a
            if ids(n) ~= 1
                tempHV = R2ind([selectedPop; PopObj(n, :)],W,r);
                if tempHV > currentHV
                    elemId = n;
                    currentHV = tempHV;
                end
            end
        end
        if elemId ~= 0
            ids(elemId) = 1;
            selectedPop = [selectedPop; PopObj(elemId, :)];
            currentHV = R2ind(selectedPop, W, r);
        end
    end


    %% output 
    Subset = selectedPop;
    res_hv = HV(Subset, r);
    res_r2 = R2ind(Subset, W, r);
end
