% 生成测试数据
clear();
rng("default");
n = 1000; m = 50; N = 20;
% n = 21000; m = 80; N = 30;
% [A, ~] = qr(rand(n, n));
% A = A(1:m,:);
% format long
A = randn(m,n);
x_true = zeros(n,1);
% x_true(1:10) = randn(10,1); % 稀疏解
% 

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

b = A*x_true + 0*randn(m,1); % 添加噪声
v = 1/(2*norm(A)^2);        % 固定步长
% norm(A*x_opt-b)
% 
% ans =
% 
%    2.0076e-07
% 
% norm(x_opt-x_true)
% 
% ans =
% 
%    2.5338e-07


% load('sglasso.mat');
% [m, n] = size(A);
% v = 1/(2*normest(A)^2);
% b = A*x;
% x_true= x;

% 算法参数
lambda_0 = 1; tau_0 = 0.1;    % 初始正则化参数
% lambda_0 = 2*1; tau_0 = 2*0.1;    % 初始正则化参数


lambda_final = 1e-4; tau_final = 1e-5; % 最终阈值
kappa = 0.96;    % 递减系数
max_iter = 1000; % 安全限制

% 运行算法
x_k = zeros(n,1); % 初始点 x^0 = 0
% x_k(1:m) = A(1:m, 1:m) \ b;
% load('xsp.mat');
% x_k = xsp;
% x_k = x_true;
[x_opt, x_history, param_history] = imtc20_n(...
    A, b, N, lambda_0, tau_0, lambda_final, tau_final, kappa, v, max_iter, x_k);
%[x_optp, x_historyp, param_historyp] = imtc2_p_n(...
%    A, b, lambda_0, tau_0, lambda_final, tau_final, kappa, v, max_iter,p);


err1 = x_true-x_opt;
%err2 = x_true-x_optp;
% norm(err1)
%norm(err2)
plot(err1,'b');
if exist('x_true')
    disp('A*x-b:');
    disp(norm(A*x_opt-b));
    disp('diff x:');
    disp(norm(x_opt-x_true));
end
xsp = x_opt;
save('xsp.mat', 'xsp');
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