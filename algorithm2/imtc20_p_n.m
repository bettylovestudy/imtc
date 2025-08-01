function [x_opt, x_history, param_history] = imtc20_p_n(A, b, lambda_0, tau_0, lambda_final, tau_final, kappa, v, max_iter,p)
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
    N = round(n / 300); % 分成20组
    group_size = floor(n / N); % 每组基本大小,向下取整
    remainder = mod(n, N); % 余数
    
    % 构建分组索引
    group_indices = cell(N, 1);
    start_idx = 1;
    for g = 1:N
        % 分配组大小，前remainder组多一个元素
        current_size = group_size + (g <= remainder);
        end_idx = start_idx + current_size - 1;
        group_indices{g} = start_idx:end_idx;
        start_idx = end_idx + 1;
    end
    
    x_history = zeros(n, max_iter+1);
    param_history = zeros(2, max_iter+1);
    
    x_k = zeros(n,1); % 初始点 x^0 = 0
    lambda_k = lambda_0;
    tau_k = tau_0;
    
    x_history(:,1) = x_k;
    param_history(:,1) = [lambda_k; tau_k];
    xdiff = 1e3;
    for k = 1:max_iter
        % 检查终止条件
        if xdiff < 1e-8 && lambda_k <= lambda_final && tau_k <=tau_final
            fprintf('算法在迭代 %d 终止：参数达到阈值 (λ=%.4f, τ=%.4f)\n', ...
                    k-1, lambda_k, tau_k);
            break;
        end
        
        % 执行迭代步骤 (3.1)-(3.3)
        % 步骤3.1: 计算y^k
        y_k = x_k - 2*v * A'*(A*x_k - b);
        
        % 步骤3.2: 计算z^k
        if p==1
            alpha = 1;
        else
            alpha = (2-p)*(2-2*p)^((p-1)/(2-p));
        end
        z_k = H_operator_p(y_k, alpha*(v*lambda_k)^(1/(2-p)),lambda_k, v, p, 1e-6, max_iter);

        
        % 步骤3.3: 计算x^{k+1} (分组处理)
        x_k_plus1 = zeros(n,1);
        for g = 1:N
            % 获取当前组的索引
            G_i = group_indices{g};
            
            % 提取z^k的对应分组
            z_Gi = z_k(G_i);
            
            % 计算分组阈值
            beta = sqrt(2*v*(lambda_k + tau_k*(norm(z_Gi,p)).^p));
            
            % 应用分组硬阈值算子
            x_k_plus1(G_i) = H_group_operator(z_Gi, beta);
        end
        
        % 更新参数和状态
        x_k = x_k_plus1;
        if lambda_k <= lambda_final || tau_k <= tau_final
            lambda_k = lambda_final;
            tau_k = tau_final;
        else
            lambda_k = kappa * lambda_k;
            tau_k = kappa * tau_k;
        end
        
        % 保存历史记录
        x_history(:,k+1) = x_k;
        param_history(:,k+1) = [lambda_k; tau_k];
        xdiff = norm(x_k-x_history(:, k));
        if mod(k, 20) == 0
            disp([k, beta, xdiff]);
        end
    end
    
    x_opt = x_k;
    
    % 截断未使用的历史记录
    x_history = x_history(:,1:k+1);
    param_history = param_history(:,1:k+1);
end