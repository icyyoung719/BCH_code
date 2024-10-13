% 由本原多项式和m求M，以构建GF(2^m)
function M = PrimitivePolynomialGenerator(m, primPoly)
    % 将十进制的本原多项式转化为二进制
    primPolyBin = de2bi(primPoly, m + 1, 'left-msb'); 

    % 将alpha初始化为一个向量，最高次幂的系数设为1
    alpha = [zeros(1, m), 1]; 

    % 初始化M数组以存储GF(2^m)的元素
    M = zeros(1, 2^m); 

    % 域中的第一个元素总是0
    M(1) = 0;

    % 使用alpha的幂生成所有非零元素
    for i = 1:(2^m - 1)
        % 将当前的alpha^i存储为域元素
        M(i + 1) = bi2de(alpha, 'left-msb');  % 将二进制向量转换为十进制数

        % 将alpha乘以x，即向左移一位
        alpha = [alpha(2:end), 0];  % 左移多项式

        % 如果最左侧的位为1，则通过本原多项式进行约减
        if alpha(1) == 1
            % 将移位后的alpha与本原多项式异或（模约减）
            alpha = xor(alpha, primPolyBin);  % 使用完整的本原多项式
        end
    end
end
