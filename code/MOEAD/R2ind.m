function R2val = R2ind(data,V,ref)
%R2 for hv approximation     
%
% data -- solution set, size n*m where n is the number of solutions and m
% is the number of objectives.
% V -- direction vectors, where each direction vector satisfies ||v||=1.
% ref -- reference point, e.g. ref = 1.1.
%
% Referene: A new R2 indicator for better hypervolume approximation. GECCO
% 2018.

    [row,dim] = size(V);
    y = 0;
    for j=1:row
        temp = abs(data-ref)./V(j,:);
        [x,~] = max(min(temp,[],2));
        y = y+x^dim;
    end
    R2val = y/row;
end