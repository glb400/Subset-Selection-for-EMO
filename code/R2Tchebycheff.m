function [R2val, mingtch] = R2Tchebycheff(data,V,ref)
%R2 with weighted Tchebycheff function    
%
% data -- solution set, size n*m where n is the number of solutions and m
% is the number of objectives.
% V -- direction vectors, where each direction vector satisfies ||v||_{1}=1.
% ref -- reference point, e.g. ref = 1.1.
%

    % [row,dim] = size(V);
    row = size(V,1);
    y = 0;
    mingtch = [];
    for j=1:row
        % temp = abs(data-ref)./V(j,:);
        temp = abs(data-ref).*V(j,:);
        [x,~] = min(max(temp,[],2));
        % y = y+x^dim;
        y = y+x;
        mingtch = [mingtch; x];
    end
    R2val = y/row;
end