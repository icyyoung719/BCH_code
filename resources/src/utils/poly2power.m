% poly2power 函数接收一个有限域 (GF) 矩阵，该多项式对应的幂次形式
function [ power ] = poly2power( GF, poly )
	% 确保多项式是 m 位宽的二进制数。
    global m;
	temp = bi2de(poly, 'left-msb');   % 将多项式从二进制转换为十进制
	poly = de2bi(temp, m, 'left-msb'); % 再将其转换为长度为 m 的二进制形式

        % 确保 poly 的列数与 GF 一致
    if size(poly, 2) < size(GF, 2)
        % 添加缺少的零，使 poly 的长度与 GF 的列数一致
        poly = [zeros(1, size(GF, 2) - size(poly, 2)), poly];
    end

	% 在 GF 中找到多项式的位置。
	location = 0;
	for i = 1:size(GF,1)
		if GF(i,:) == poly  % 如果在 GF 中找到与 poly 相同的元素
			location = i;    % 记录位置
		end
	end

	power = location - 2;  % 幂次形式的值是位置减去 2
end
