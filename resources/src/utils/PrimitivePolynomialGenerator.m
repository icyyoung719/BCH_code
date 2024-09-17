% 由本原多项式和m求M，以构建GF(2^m)
function M = PrimitivePolynomialGenerator(m)
    % primPolys存放不同m对应的本原多项式的10进制表示 
    primPolys = [
        7,      % 对应m=2，X^2 + X + 1 -> 111 (binary) = 7 (decimal)
        11,     % 对应m=3，X^3 + X + 1 -> 1011 (binary) = 11 (decimal)
        19,     % 对应m=4，X^4 + X + 1 -> 10011 (binary) = 19 (decimal)
        37,     % 对应m=5，X^5 + X^2 + 1 -> 100101 (binary) = 37 (decimal)
        67,     % 对应m=6，X^6 + X + 1 -> 1000011 (binary) = 67 (decimal)
        131,    % 对应m=7，X^7 + X + 1 -> 10000011 (binary) = 131 (decimal)
        285,    % 对应m=8，X^8 + X^4 + X^3 + X^2 + 1 -> 100011101 (binary) = 285 (decimal)
        529,    % 对应m=9，X^9 + X^4 + 1 -> 1000010001 (binary) = 529 (decimal)
        1033,   % 对应m=10，X^10 + X^3 + 1 -> 10000001001 (binary) = 1033 (decimal)
        2053,   % 对应m=11，X^11 + X^2 + 1 -> 100000000101 (binary) = 2053 (decimal)
        4179    % 对应m=12，X^12 + X^6 + X^4 + X + 1 -> 100000100011 (binary) = 4179 (decimal)
    ];
    
    % 由m选择预先准备的本原多项式
    primPoly = primPolys(m - 1);  % 根据m选择相应的本原多项式
    % Convert the primitive polynomial from decimal to binary
    primPolyBin = de2bi(primPoly, m + 1, 'left-msb');  % Get the binary form of the polynomial

    % Initialize the alpha value as a vector with the highest power coefficient set to 1
    alpha = [zeros(1, m) , 1];  % Corresponds to alpha (the generator)

    % Initialize the M array to store elements of GF(2^m)
    M = zeros(1, 2^m);  % There are 2^m elements in GF(2^m)

    % First element in the field is always 0
    M(1) = 0;

    % Generate all non-zero elements using powers of alpha
    for i = 1:(2^m - 1)
        % Store the current alpha^i as a field element
        M(i+1) = bi2de(alpha, 'left-msb');  % Convert binary vector to decimal number

        % Multiply alpha by x, i.e., left shift by one position
        alpha = [alpha(2:end),0];  % Left shift the polynomial

        % If the leftmost bit is 1, reduce by primitive polynomial
        if alpha(1) == 1
            % XOR the shifted alpha with the primitive polynomial (modular reduction)
            alpha = xor(alpha, primPolyBin);  % Use the full primitive polynomial
        end
    end
end