function result = Ncut_WithoutSimilarityMatrix(P, H, label)

% H, P: N*k
% label: N*1
% 复杂度: O(kN^2)
    
    [n, ~] = size(H);
    classNum = length(unique(label));
    result = 0;
    
    for i = 1:classNum
        
        Vl = find(label == i);
        m = length(Vl);
        
        r = zeros(m, n);
        for j = 1:m
            % P(Vl(j), :): 1*k
            W_ = P(Vl(j), :) * H'; % 复杂度: O(Nk)
            r(j, :) = W_;
        end
        V = sum(sum(r));
        r(:, Vl) = 0;
        s = sum(sum(r));
        result = result + s / V;
        
    end 

end