function y = mix_threshold(x, t, alpha)
% MIX_THRESHOLD 混合阈值函数，结合L0和L1阈值特性
% 输入:
%   x - 输入向量
%   t - 阈值参数
%   alpha - 混合参数
% 输出:
%   y - 阈值处理后的向量

    % L1阈值 (软阈值)
    soft_threshold = sign(x) .* max(abs(x) - (1-alpha)*t, 0);
    
    % L0阈值 (硬阈值)
    hard_threshold = x .* (abs(x) > alpha*t);
    
    % 混合两种阈值结果
    y = soft_threshold .* (abs(x) > alpha*t) + hard_threshold .* (abs(x) <= alpha*t);
end

