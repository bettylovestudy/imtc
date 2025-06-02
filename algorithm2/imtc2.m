function [x_opt, x_all, param_all] = imtc2(A, b, lambda_0, tau_0, lambda_final, tau_final, kappa, v, max_iter)
    % 输入参数：
    % A: 矩阵 A
    % b: 观测向量 b
    % lambda_0, tau_0: 初始正则化参数
    % lambda_final, tau_final: 最终正则化参数阈值
    % kappa: 递减系数 (0 < kappa < 1)
    % v: 固定步长
    % max_iter: 最大迭代次数（防止无限循环）
    
    % 初始化
    n = size(A, 2);
    N = n; % 假设每组只有一个元素，实际应用中根据分组结构调整
    x_all = zeros(n, max_iter+1);%保存所有迭代过程中的x
    param_all = zeros(2, max_iter+1);
    
    x_k = zeros(n,1); % 初始点 x^0 = 0
    lambda_k = lambda_0;
    tau_k = tau_0;
    
    x_all(:,1) = x_k;
    param_all(:,1) = [lambda_k; tau_k];
    
    for k = 1:max_iter
        % 检查终止条件
        if lambda_k < lambda_final || tau_k < tau_final
            fprintf('算法在迭代 %d 终止：参数达到阈值 (λ=%.4f, τ=%.4f)\n', ...
                    k-1, lambda_k, tau_k);
            break;
        end
        
        % 执行迭代步骤 (3.1)-(3.3)
        % 步骤3.1: 计算y^k
        y_k = x_k - 2*v * A'*(A*x_k - b);
        
        % 步骤3.2: 计算z^k
        z_k = H_operator(y_k, sqrt(2*v*tau_k));
        
        % 步骤3.3: 计算x^{k+1} (分组处理)
        x_k_1 = zeros(n,1);
        for i = 1:N
            G_i = i; % 简化处理：假设每组只有一个元素
            z_Gi = z_k(G_i);
            
            % 计算分组阈值
            beta = sqrt(2*v*(lambda_k + tau_k*nnz(z_Gi)));
            
            % 应用分组硬阈值算子
            x_k_1(G_i) = H_group_operator(z_Gi, beta);
        end
        
        % 更新参数和状态
        x_k = x_k_1;
        lambda_k = kappa * lambda_k;
        tau_k = kappa * tau_k;
        
        % 保存历史记录
        x_all(:,k+1) = x_k;
        param_all(:,k+1) = [lambda_k; tau_k];
    end
    
    x_opt = x_k;
    
    % 截断未使用的历史记录
    x_all = x_all(:,1:k+1);
    param_all = param_all(:,1:k+1);
end