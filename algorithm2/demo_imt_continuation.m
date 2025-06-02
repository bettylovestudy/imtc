
% 设置参数
n = 1000;       % 信号维度
m = 200;        % 测量数
k = 50;         % 稀疏度
alpha = 0.3;    % 混合参数
lambda_max = 1.0;
lambda_min = 0.01;
iter_max = 200;
tol = 1e-6;

% 生成稀疏信号
x_true = zeros(n, 1);
idx = randperm(n, k);
x_true(idx) = 5 * randn(k, 1);

% 生成测量矩阵
A = randn(m, n);
A = A ./ sqrt(sum(A.^2, 1));  % 归一化列向量

% 生成观测向量
b = A * x_true + 0.1 * randn(m, 1);  % 添加噪声

% 运行IMT算法
[x, history] = imt_continuation(A, b, alpha, lambda_max, lambda_min, iter_max, tol);

% 评估结果
recovery_error = norm(x - x_true) / norm(x_true);
fprintf('信号恢复误差: %.4f\n', recovery_error);
fprintf('估计的稀疏度: %d\n', nnz(x));
fprintf('真实的稀疏度: %d\n', k);

% 可视化结果
figure;
subplot(2, 1, 1);
stem(1:n, x_true, 'b', 'MarkerSize', 3);
title('真实信号');
subplot(2, 1, 2);
stem(1:n, x, 'r', 'MarkerSize', 3);
title('恢复信号');

figure;
semilogy(history.iter, history.obj, 'b-', 'LineWidth', 1.5);
xlabel('迭代次数');
ylabel('目标函数值');
title('收敛曲线');
grid on;
