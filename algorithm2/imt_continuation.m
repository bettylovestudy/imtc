function [x, history] = imt_continuation(A, b, alpha, lambda_max, lambda_min, iter_max, tol)
% IMT_CONTINUATION 带有连续化技术的迭代混合阈值算法
% 用于求解混合稀疏优化问题: min_x 0.5*||Ax-b||^2 + lambda*(alpha*||x||0 + (1-alpha)*||x||1)
% 输入:
%   A - 测量矩阵 (m×n)
%   b - 观测向量 (m×1)
%   alpha - 混合参数 (0≤alpha≤1)，控制L0和L1的权重
%   lambda_max - 初始正则化参数
%   lambda_min - 最终正则化参数
%   iter_max - 最大迭代次数
%   tol - 收敛容差
% 输出:
%   x - 稀疏解 (n×1)
%   history - 迭代历史记录

    % 初始化
    [m, n] = size(A);
    x = zeros(n, 1);
    lambda = lambda_max;
    mu = 1.0 / norm(A'*A, 2);  % 步长参数
    decay_factor = (lambda_min / lambda_max)^(1/iter_max);  % 正则化参数衰减因子
    
    % 记录历史
    history.iter = [];
    history.obj = [];
    history.residual = [];
    history.nnz = [];
    
    % 主循环
    for iter = 1:iter_max
        % 梯度步骤
        grad = A' * (A * x - b);
        x_temp = x - mu * grad;
        
        % 混合阈值操作
        x_new = mix_threshold(x_temp, lambda * mu, alpha);
        
        % 更新正则化参数
        lambda = lambda * decay_factor;
        
        % 计算目标函数值
        obj = 0.5 * norm(A * x_new - b)^2 + lambda * (alpha * nnz(x_new) + (1-alpha) * norm(x_new, 1));
        
        % 记录历史
        history.iter = [history.iter; iter];
        history.obj = [history.obj; obj];
        history.residual = [history.residual; norm(A * x_new - b)];
        history.nnz = [history.nnz; nnz(x_new)];
        
        % 检查收敛性
        if norm(x_new - x) / norm(x) < tol
            break;
        end
        
        x = x_new;
    end
end

