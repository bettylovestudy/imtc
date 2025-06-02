n = 1000; m = 50;
[A, ~] = qr(rand(n, n));
A = A(1:m,:);
% A = randn(m,n);
x_true = zeros(n,1);
% x_true(1:10) = randn(10,1); % 稀疏解
N = 20;
gn = floor(n / N);
s = 0;
S = 0;
for i = 1 : N
    shift = (i - 1) * gn;
    num = rand;
    if i == 1 || i == N
        num = 0.96;
    end
    if num > 0.95
        S = S + 1;
        disp([i, shift, num]);
        for j = 1 : (gn / 10)
            ind = randi(gn);
            x_true(shift+ind) = rand;
            s = s + 1;
        end
    end
    % disp(x_true(shift+1:shift+gn)');
end
A = A / (norm(A*x_true) / norm(x_true));

disp(['s: ', num2str(s)]);
disp(['S: ', num2str(S)]);
norm_Ax = norm(A*x_true);
norm_x = norm(x_true);
delta = max(1-norm_Ax^2/norm_x^2, norm_Ax^2/norm_x^2-1);
disp(['delta: ', num2str(delta)]);

lambda_0 = 1; tau_0 = 0.1;    % 初始正则化参数
lambda_final = 1e-4; tau_final = 1e-5; % 最终阈值
kappa = 0.96;    % 递减系数
v = 1/(2*norm(A)^2);        % 固定步长
flag = delta*(1+sqrt(s)+sqrt(S*(lambda_0/tau_0+ceil(n/N))));
eta = (1 - flag) / 2;