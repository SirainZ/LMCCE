function out = mean_WithoutSimilarityMatrix(x, num, sigma)

[n, ~] = size(x);
finalnum = mod(n, num);

if finalnum ~= 0
    m = (n - finalnum) / num + 1;
    
    Y = zeros(n, 1);
    for t = 1:m-1       
        xt = x(num*t-num+1:num*t, :);
        
        D = pdist2(xt, x);
        Rt = exp(-(D / sigma).^2);
        
        Yt = mean(Rt, 2);
        Y(num*t-num+1:num*t, :) = Yt;
    end
    
    xt = x(n-finalnum+1:n, :);
    
    D = pdist2(xt, x);
    Rt = exp(-(D / sigma).^2);
    
    Yt = mean(Rt, 2);
    Y(n-finalnum+1:n, :) = Yt;
    
    
else
    m = n / num;
    
    Y = zeros(n, 1);
    for t = 1:m       
        xt = x(num*t-num+1:num*t, :);
        
        D = pdist2(xt, x);
        Rt = exp(-(D / sigma).^2);
              
        Yt = mean(Rt, 2);
        Y(num*t-num+1:num*t, :) = Yt;
    end
    
end

out = mean(Y);

end