function A = Orth_Matrix(m,n)
A = randn(m,n);
for i = 1:m
    % 当前行
    v = A(i, :);
    % 减去与前面所有正交行的投影
    for j = 1:i-1
        v = v - (A(j, :) * v') * A(j, :);
    end
    % 归一化，使行向量的模长为1
    norm_v = norm(v);
    if norm_v > 1e-10  % 避免除以零
        A(i, :) = v / norm_v;
    else
        error('矩阵行向量线性相关，无法正交化');
    end
end
end