function v = LMCCE_cnum(H, kmin, interval, kmax)

    [n, ~] = size(H);

    klist = kmin: interval: kmax;
    cnum = zeros([length(klist), 1]);
    cnum(1) = n;
    cnum(2) = n;

    S = H' * H;
    Sk = S ^ klist(1);
    Ss = S ^ interval;
    for i = 3:length(klist)
        Sk = Sk * Ss;  % S^3=HS^2H'
    
        idx = zeros(n, 1);
        P = H * Sk;
        for j = 1:n
            Wi = P(j, :) * H';
            [~, idx(j)] = max(Wi, [], 2);
        end
        c = find((1:n) == idx')';
        cnum(i) = length(c);   
    end

    v = [klist', cnum];

end