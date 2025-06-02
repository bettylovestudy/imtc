% 硬阈值算子H 
function z = H_operator(y, alpha)
    z = zeros(size(y));
    y_star = abs(y) > alpha;
    z(y_star) = y(y_star);
end

