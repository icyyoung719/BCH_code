% 定义有限域元素
M = [0 1 2 4 8 16 5 10 20 13 26 17 7 14 28 29 31 27 19 3 6 12 24 21 15 30 25 23 11 22 9 18];
GF = de2bi(M, 5, 'left-msb');

% 计算最小多项式
function minimal_poly = compute_minimal_poly(GF, alpha)
    n = length(alpha);
    minimal_poly = poly(ones(n, 1));
    
    for i = 1:n
        x_i = GF(i) + 1;
        
        % 使用辗转相除法求解最小多项式
        while rem(minimal_poly, x_i) == 0
            minimal_poly = x_i * minimal_poly / gcd(minimal_poly, x_i);
        end
        
        if minimal_poly == 1
            break;
        end
    end
    
    minimal_poly = poly(minimal_poly);
end

% 示例：计算第一个元素的最小多项式
alpha = GF(1,:);
minimal_poly_1 = compute_minimal_poly(GF, alpha)

% 输出结果：
% minimal_poly_1 =
%     1     1     0