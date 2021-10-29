function [x_new] = update_my(x, alpha, delta_x)
    for i = 1:1:length(x)
        for j = 1:1:length(x{i})
            x_new{i}{j} = x{i}{j} - alpha*delta_x{i}{j};
        end
    end
end

