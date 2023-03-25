function [map, y] = label_corrected(x, gt)

    % gt, 按各类地物像素点数排序
    n = length(unique(gt));
    t = tabulate(gt);
    t = sortrows(t, 2,'descend');
    t = t(1:n, 1);

    % 按块对标签
    x = x + 100;  
    u = unique(x);
    map = zeros(2, n);
    for i = 1:n
        p = find(gt == t(i));
        part = x(p);
        temp = sortrows(tabulate(part), 2,'descend');
        temp = temp(:, 1);
        rem = ismember(temp, map(1, :));
        sel = find(rem == 0);
        sel = temp(sel(1));
        
        map(1, i) = sel;
        map(2, i) = t(i);                 
    end
    
    % 得到map
    [~, r] = find(map(1, :) < 100);
    los = setdiff(u, map(1, :));
    s = zeros(2, length(los));
    s(1, :) = los;
    for i = 1:length(los)
        s(2, i) = length(find(x == los(i)));
    end
    s = sortrows(s', 2, 'descend');
    s = s(:, 1);
    
    for i = 1:length(los)
        map(1, r(i)) = s(i);
    end
    
    
    % 修正label
    for i = 1:n
        x(x == map(1, i)) = map(2, i);
    end
    y = x;
    
end