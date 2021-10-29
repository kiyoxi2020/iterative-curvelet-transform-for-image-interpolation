
function y = threshold_my(x, lambda, type)
if nargin < 3
    type = 'soft';
end
switch type
    case 'soft'
        for i = 1:1:length(x)
            for j = 1:1:length(x{i})
                x{i}{j}   = (max(abs(x{i}{j}) - lambda, 0)).*sign(x{i}{j});
            end
        end
        
    case 'hard'
        for i = 1:1:length(x)
            for j = 1:1:length(x{i})
                x{i}{j}   = (abs(x{i}{j}) > lambda).*x{i}{j};
            end
        end
end

y = x;
end