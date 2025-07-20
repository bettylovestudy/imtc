% 生成测试数据
rng("default");
n = 1000; m = 50;
[A, ~] = qr(rand(n, n));
A = A(1:m,:);
% A = randn(m,n);
x_true = zeros(n,1);
% x_true(1:10) = randn(10,1); % 稀疏解
N = 20;
gn = floor(n / N);
s = 0;
for i = 1 : N
    shift = (i - 1) * gn;
    num = rand;
    if i == 1 
        num = 0.96;
    end
    if num > 0.95
        disp([i, shift, num]);
        for j = 1 : 2
            ind = randi(gn);
            x_true(shift+ind) = rand;
            s = s + 1;
        end
    end
    % disp(x_true(shift+1:shift+gn)');
end
A = A / (norm(A*x_true) / norm(x_true));
% disp([norm(x_true), norm(A*x_true)]);
disp(['nonzero ratio: ', num2str(s/n)]);

b = A*x_true + 1e-3*randn(m,1); % 添加噪声
%p=1/2;%给p赋值
% 算法参数
lambda_0 = 1; tau_0 = 0.1;    % 初始正则化参数

lambda_final = 1e-4; tau_final = 1e-5; % 最终阈值
kappa = 0.96;    % 递减系数
v = 1/(2*norm(A)^2);        % 固定步长
max_iter = 1000; % 安全限制

% 运行算法
x_k = zeros(n,1); % 初始点 x^0 = 0
% x_k = x_true;
[x_opt2p, x_historyp, param_historyp] = imtc22_p_n(...
   A, b, lambda_0, tau_0, lambda_final, tau_final, kappa, v, max_iter,p);

err1 = x_true-x_opt2p;
%err2 = x_true-x_optp;
norm(err1)
%norm(err2)
plot(err1,'b');
%hold on
%plot(err2(1:30),'r')

% f = @(x) 1/2*norm(A*x-b)+lambda_final*sum() +tau_final*sum(abs(x)^0);
% 绘制结果
% figure;
% subplot(2,1,1);
% plot(x_true, 'ro', 'DisplayName', '真实解'); hold on;
% plot(x_opt, 'b*', 'DisplayName', '算法解'); 
% legend(); title('解的比较');
% 
% subplot(2,1,2);
% semilogy(param_history(1,:), 'r-', 'DisplayName', 'λ_k'); hold on;
% semilogy(param_history(2,:), 'b--', 'DisplayName', 'τ_k');
% legend(); title('参数递减过程');
% xlabel('迭代次数');