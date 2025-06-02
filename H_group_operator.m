% 定义分组硬阈值算子H_{G_i} (3.5)
function x = H_group_operator(z, beta)
    if norm(z) > beta
        x = z;
    else
        x = zeros(size(z));
    end
end