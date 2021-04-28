function [selected] = Closest_selection(selected, N, D, scores)
%CLOSEST_SELECTION Summary of this function goes here
%   Detailed explanation goes here

next = ones(size(scores, 1), 1);

[sD, sDI] = sort(D, 1);

P = size(selected,2);

k = sum(selected) - N;

while k > 0
    min_i = 0;
    min_j = 0;
    min_val = inf;
    
    for i=1:P
        if selected(i)
            next_i = next(i);
            ii = sDI(next_i, i);
            while ~selected(ii)
                next(i) = next_i + 1;
                next_i = next(i);
                ii = sDI(next_i, i);
            end
            i_val = sD(next_i, i);
            if i_val < min_val
                min_i = i;
                min_j = ii;
                min_val = i_val;
            end
        end
    end
    
    score_i = scores(min_i);
    score_j = scores(min_j);
    
    to_eliminate = min_i;
    
    if score_i < score_j
        to_eliminate = min_j;
    end
    
    selected(to_eliminate) = false;
    
    k = k - 1;
end


end

