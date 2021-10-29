function [y] = norm_my_l1(x)
    y = 0;
    for i = 1:1:length(x)
        for j = 1:1:length(x{i})
            y = y + sum(sum(abs(x{i}{j})));
        end
    end
end

