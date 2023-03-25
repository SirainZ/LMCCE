function [OA, Kappa, PA, UA] = calaccuracy(gt, ind)

    n = length(ind);

    mat = confusionmat(gt, ind);
    c = mat';
    
    % OA
    p0 = sum(diag(mat)) / n;
    OA = p0;
    
    % Kappa
    fact = sum(c, 1);
    pred = sum(c, 2);
    pe = sum(fact .* pred') / (n^2);
    Kappa = (p0 - pe) / (1 - pe);
    
    % PA
    k = unique(gt);
    PA = zeros(length(k), 1);
    for i = 1:length(k)
        PA(i) = c(i, i) / sum(c(:, i));
    end
    
    % UA
    UA = zeros(length(k), 1);
    for i = 1:length(k)
        UA(i) = c(i, i) / sum(c(i, :));
    end

end