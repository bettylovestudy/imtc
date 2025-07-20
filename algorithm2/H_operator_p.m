% 硬阈值算子H_p 
function z = H_operator_p(y, alpha, lambda, v, p, tol, max_iter)
    n = length(y);
    z = zeros(n, 1);
    for i = 1:n
        yi_abs = abs(y(i));
        if yi_abs <= alpha
            z(i) = 0;
        else
            % |z| > h_a > beta_a
            % abs(y(i)) > alpha > a_lower
            assert(yi_abs >= alpha);
            
            % 初始值选择：a^0 ∈ [(2λ_k v_k (1-p))^(1/(2-p)), |y_i|]
            a_lower = (2 * lambda * v * (1 - p))^(1/(2 - p));
            assert(alpha >= a_lower);
            a0 = min(max(a_lower, 0.5 * yi_abs), yi_abs);  % 取中间值作为初始猜测
            
            % 迭代求解 a*
            a_prev = a0;
            converged = false;
            
            for iter = 1:max_iter
                % 迭代公式
                a_next = yi_abs - lambda * v * p * a_prev^(p - 1);
                assert(a_next>=a_prev);
                assert(a_next<=yi_abs);
                % 检查收敛性
                if abs(a_next - a_prev) < tol
                    converged = true;
                    break;
                end
                
                a_prev = a_next;
            end
            if ~converged
                warning(['迭代达到最大次数 ', num2str(max_iter), ...
                       ' 但未收敛到容差 ', num2str(tol), '。']);
            end
            % 使用收敛值或最后一次迭代结果
            a_star = a_prev;
            z(i) = sign(y(i)) * a_star; 
        end
    end
end

