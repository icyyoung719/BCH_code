function root_sigma_X = find_root_sigma_X(sigma_X,m,field_table)
    % find_root_sigma_X 函数用于计算错误定位多项式的根，即错误的位置
    % - sigma_X: 错误定位多项式
    % - m: 有限域 GF(2^m) 的次幂
    % - field_table: 有限域 GF(2^m) 的查找表，用于执行加法和乘法操作
    % - root_sigma_X: 错误定位多项式的根（即错误的位置）

	root_sigma_X = [];
    % 遍历有限域中的所有元素，从 α^1 到 α^(2^m-1)，依次代入检验
	for i = 2:2^m
		tmp = [0];
        % 计算错误定位多项式在当前值 (i-2) 处的值
		for j = 2:length(sigma_X)
			if sigma_X(j) ~= -1
                % 表示 sigma_X 多项式的每一项的值
				tmp(end+1) = mod(sigma_X(j)+(i-2)*(j-1),2^m-1);
			end
		end

        % 对多项式的每一项逐步进行有限域加法，得到多项式tmp的值
		for j = 1:length(tmp)-1
			tmp(j+1) = gf_add(tmp(j),tmp(j+1),field_table);
		end

        % 如果多项式tmp值为 -1（即等于0），则 i-2 是多项式的根
		if tmp(end) == -1
			root_sigma_X(end+1) = i-2;
		end
	end
end