function [sigma_x] = berlekamp_massey(S_a, t, GF, primPolyBin)
    % 输入：
    % S_a: 接收到的码字的综合 (syndrome) 序列，初始化类型为 S_a = zeros(1, 2*t);
    % t: 纠错能力，最大纠错的错误个数
    % GF: 有限域的二进制表示 GF = de2bi(M, m, 'left-msb');
    % primPolyBin: primPolyBin = de2bi(primPoly, m + 1, 'left-msb');  % 本原多项式的二进制形式
    % Initialize variables
    global m;
    global k;
    n = 2^m - 1;  % Length of the codeword
    mu = zeros(1, 2*t);
    mu(1) = -1/2;  % Initial value of mu

    % Sigma matrix initialization
    sigma_x = zeros(n, n, 2*t);  % Assuming n is the max degree we might deal with
    sigma_x(1, 1, 1) = 1;  % Start with sigma_0 = 1

    d_mu = zeros(1, 2*t);
    l_mu = zeros(1, 2*t);

    twomu_lmu = zeros(1, 2*t);
    twomu_lmu(1) = -1;  % Initial value

    % Fill in the table
    for j = 1:2*t
        if j > 1
            d_mu(j) = S_a(j-1);  % Fill d_mu with syndrome values
        end

        if d_mu(j) == 0
            % No update for sigma if d_mu is zero
            sigma_x(:,:,j+1) = sigma_x(:,:,j);
        else
            % Find a preceding row with the most positive twomu_lmu and d_mu ~= 0
            rho = 0;
            most_pos = -Inf;

            for k = 1:j-1
                if twomu_lmu(k) > most_pos && d_mu(k) ~= -1
                    most_pos = twomu_lmu(k);
                    rho = k;
                end
            end

            % Calculate the new sigma_x
            power = 2 * (mu(j) - mu(rho));
            sigma_rho = sigma_x(:,:,rho);

            % Shift sigma_rho according to the power
            for i = 1:power
                sigma_rho(1:end-1) = sigma_rho(2:end);
                sigma_rho(end) = 0;  % Assuming -1 is treated as 0 in binary field
            end

            coeff = p_mult(d_mu(j), p_inv(d_mu(rho)));  % Calculate d_mu * d_rho^(-1)

            % Update sigma_rho with the calculated coefficient
            for i = 1:length(sigma_rho)
                if sigma_rho(i) ~= -1
                    sigma_rho(i) = p_mult(sigma_rho(i), coeff);
                end
            end

            % Update sigma_x with the new values
            sigma_x(:,:,j+1) = poly_add(GF, sigma_x(:,:,j), sigma_rho);
        end

        % Terminate if mu reaches error correction limit
        if mu(j) == t - 1
            break;
        end

        % Calculate l_mu(j+1)
        degree = max(find(sigma_x(:,1,j+1) ~= -1));
        l_mu(j+1) = degree;

        % Update d_mu(j+1) using S_a and sigma_x
        S_sigma_coeff = zeros(2, 2);
        for i = 0:l_mu(j+1)
            S_sub = 2 * mu(j) + (3 - i);
            sigma_term = length(sigma_x(:,:,j+1)) - i;
            S_sigma_coeff(:, i + 1) = [S_sub; sigma_term];
        end

        alpha_poly = -1 * ones(1, 3);
        for i = 1:length(S_sigma_coeff)
            if S_sigma_coeff(1, i) ~= -1
                alpha_poly(i) = p_mult(S_a(S_sigma_coeff(1, i)), sigma_x(1, S_sigma_coeff(2, i), j+1));
            end
        end

        % Create binary representation of alpha polynomial
        binary_alpha_poly = zeros(1, 500);
        for i = 1:length(alpha_poly)
            if alpha_poly(i) ~= -1
                binary_alpha_poly(length(binary_alpha_poly) - alpha_poly(i)) = 1;
            end
        end

        % Update d_mu for the next iteration
        d_mu(j+1) = poly2power(GF, polynomial_mod(binary_alpha_poly, primPolyBin));

        % Calculate twomu_lmu(j+1)
        twomu_lmu(j+1) = 2 * mu(j+1) - l_mu(j+1);
    end

    % Error locator polynomial is the last entry of sigma_x
    error_loc_poly = sigma_x(:,:,length(sigma_x));
end