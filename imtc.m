function [x_opt, x_history] = imtc(A, b, x0, lambda, tau, lambda_k, tau_k, v_k, max_iter)
    % 输入参数：
    % A: 矩阵 A
    % b: 观测向量 b
    % x0: 初始点 x^0
    % lambda, tau: 最终参数值
    % lambda_k, tau_k: 递减的参数序列（应为长度为max_iter的向量）
    % v_k: 步长序列（应为长度为max_iter的向量）
    % max_iter: 最大迭代次数
    
    % 初始化
    n = length(x0);
    N = size(A, 2); % 假设N是A的列数，对应分组数
    x_history = zeros(n, max_iter+1);
    x_history(:,1) = x0;
    x_k = x0;
    
    for k = 1:max_iter
        % 步骤3.1: 计算y^k
        y_k = x_k - 2*v_k(k) * A'*(A*x_k - b);
        
        % 步骤3.2: 计算z^k
        z_k = H_operator(y_k, sqrt(2*v_k(k)*tau_k(k)));
        
        % 步骤3.3: 计算x^{k+1} (分组处理)
        x_k_plus1 = zeros(n,1);
        for i = 1:N
            % 假设G_i是第i组索引，这里简化处理为单个元素
            % 实际应用中需要根据具体分组情况调整
            G_i = i; % 这里假设每组只有一个元素，实际应替换为真实分组
            
            % 提取z^k的对应分组
            z_Gi = z_k(G_i);
            
            % 计算分组阈值
            beta = sqrt(2*v_k(k)*(lambda_k(k) + tau_k(k)*nnz(z_Gi)));
            
            % 应用分组硬阈值算子
            x_k_plus1(G_i) = H_group_operator(z_Gi, beta);
        end
        
        % 更新x_k并保存历史
        x_k = x_k_plus1;
        x_history(:,k+1) = x_k;
    end
    
    x_opt = x_k;
end

