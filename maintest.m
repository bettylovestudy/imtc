% 生成测试数据
n = 100; m = 50;
A = randn(m,n);
x_true = zeros(n,1);
x_true(1:10) = randn(10,1); % 稀疏解
b = A*x_true + 0.1*randn(m,1); % 添加噪声

% 算法参数
x0 = zeros(n,1); % 初始点
max_iter = 100;
lambda = 0.1; tau = 0.1;
lambda_k = linspace(1, lambda, max_iter); % 递减序列
tau_k = linspace(1, tau, max_iter); % 递减序列
v_k = 0.1*ones(max_iter,1); % 固定步长

% 运行算法
[x_opt, x_history] = imtc(A, b, x0, lambda, tau, lambda_k, tau_k, v_k, max_iter);

% 绘制结果
figure;
plot(x_true, 'ro', 'DisplayName', '真实解');
hold on;
plot(x_opt, 'b*', 'DisplayName', '算法解');
legend();
title('算法解与真实解对比');