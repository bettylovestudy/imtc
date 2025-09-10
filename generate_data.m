%% 模拟数据A
clear;clc;
rng("default");
% 生成m×n的标准正态分布随机矩阵
m = 256;
n = 1024;
N = 64;
tol = 1e-16;
N_expe = 10;
rel_error = zeros(N_expe,1);
for iter = 1 : N_expe
    % 对矩阵进行正交化处理，满足A*A'=I
    A = Orth_Matrix(m,n);
    % 模拟数据x
    % 创建gLen×N的零矩阵
    gLen = round(n/N);
    assert(gLen*N == n);
    x = zeros(gLen, N);
    % 随机选择S个组
    S = 20;
    Gind = randperm(N, S);

    % 为选中的组填充随机值
    for u = 1:length(Gind)
        % 每个组中随机选择r个元素
        r = 1;
        ind = randperm(gLen, r);
        % 为选中位置赋值标准正态分布随机数
        x(ind, Gind(u)) = randn(r, 1);
        % disp(x(ind, Gind(u)));
    end

    % 将x转换为单列矩阵
    true_x = x(:);

    % 模拟数据b
    sigma = 0.001;
    % 生成噪声项
    epsilon = sigma * randn(m, 1);
    % 计算带噪声的响应变量
    b = A * true_x + epsilon;

    %% 算法参数
    v = 0.5;        % 固定步长
    lambda_0 = 1; tau_0 = 0.1;    % 初始正则化参数

    lambda_final = 1e-4*lambda_0;
    tau_final = 1e-5*tau_0;% 最终阈值
    kappa = 0.98;    % 递减系数
    max_iter = 1000; % 安全限制

    %% 运行算法
    x_k = zeros(n,1); % 初始点 x^0 = 0
    % time0 = datetime("now");
    [x_opt, x_history, param_history] = imtc20_n(...
        A, b, N, lambda_0, tau_0, lambda_final, tau_final, kappa, v, max_iter, x_k);
    % time1 = datetime("now") - time0;

    %% error
    err1 = x_opt - true_x;

    if exist('true_x')
        disp(['A*x-b:imtc ',num2str(norm(A*x_opt-b))]);
        rel_error(iter) = norm(err1)/norm(true_x);
        % disp(['rel_error x:',num2str(rel_error(iter))]);
        % disp(['sparse x: ',num2str(nnz(x_opt))]);
        % disp(['time: ',num2str(time1)]);
    end
end
disp('************************************************************');
% disp(['sparse ratio: ',num2str(nnz(true_x)/n)]);
% disp('on intra- and inter-group sparsity: ');
fprintf("组内组间稀疏比率（gamma）: %.2f%%\n", r/S * 100);
disp(['lambda_0: ', num2str(lambda_0), ', tau_0: ', num2str(tau_0)]);
disp('below 0.005 rate(imtc): ')
disp(sum(rel_error < 0.005)/N_expe);
disp('relative error details: ')
disp(rel_error');

