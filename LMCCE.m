function [label_out, label_all, sel, ncut_all, cnum] = LMCCE(H, kmin, interval, kmax, classnum)

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

loc = find(cnum == classnum);
l = length(loc);
label_all = zeros(n, l);
ncut_all = zeros(1, l);
for i = 1:length(loc)
    Sk = S ^ (klist(loc(i)-1));
    idx = zeros(n, 1);
    P = H * Sk;
    for j = 1:n
        Wi = P(j, :) * H';
        [~, idx(j)] = max(Wi, [], 2);
    end
    c = find((1:n) == idx')';
    m = length(c);
    
    Q = zeros(m, n);
    for j = 1:m
        W_ = P(c(j), :) * H';
        Q(j, :) = W_;
    end
    
    [~, label] = max(Q, [], 1);
    label = label';
    %         centers = c;
    
    P_ori = H * (S ^ 0);
    ncut = Ncut_WithoutSimilarityMatrix(P_ori, H, label);
    
    label_all(:, i) = label;
    ncut_all(i) = ncut;
end

[~, sel] = min(ncut_all);
sel = sel(1);
label_out = label_all(:, sel);

end


