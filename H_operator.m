% 定义硬阈值算子H (3.4)
function z = H_operator(y, alpha)
    z = zeros(size(y));
    mask = abs(y) > alpha;
    z(mask) = y(mask);
end

