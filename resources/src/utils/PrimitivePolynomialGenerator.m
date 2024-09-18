% 由本原多项式和m求M，以构建GF(2^m)
function M = PrimitivePolynomialGenerator(m,primPoly)
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