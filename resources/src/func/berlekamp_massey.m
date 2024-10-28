function sigma_X = berlekamp_massey(syndrome,t,m,field_table)
    % 用于BCH码的错误定位多项式计算
    % - syndrome: 伴随式向量，用于计算错误定位多项式
    % - t: 纠错能力
    % - m: GF(2^m) 域的次幂
    % - field_table: 有限域 GF(2^m) 的查找表

    % 初始化 BM_table，用于记录算法的迭代过程
    % - 第1列: 迭代步索引
    % - 第2列: 错误定位多项式的当前项表示
    % - 第3列: 差值 d(k)
    % - 第4列: 错误定位多项式的当前阶数
    % - 第5列: 2*(k+1) - 错误定位多项式的阶数
	BM_table = cell(t+2,5);				
    %the second column of RM_table is an one-dimensions array, this row is the order of alpha
	%ie. 1+X+(alpha^5)*X^2
	%[0 0 5]
	%ie. 1+X+(alpha^5)*X^3
	%[0 0 -1 5]		%-1 denote this X^i do not exist
    
    % 初始化第一行，表示 k = -1/2 的初始状态
	BM_table{1,1} = -1/2;	BM_table{1,2} = [0]; BM_table{1,3} = 0; BM_table{1,4} = 0; BM_table{1,5} = -1;
    % 初始化第二行，表示 k = 0 的状态
	BM_table{2,1} = 0;	BM_table{2,2} = [0]; BM_table{2,4} = 0; BM_table{2,5} = 0;

    % 迭代循环，从 k = 0 开始，到 k = t - 1 为止
	for k_ = 0:t-1
		dk = compute_dk(BM_table, syndrome, field_table, m, k_);
        % 如果差值 d(k) 为 -1，说明没有错误的校正项，继续沿用上一行的多项式
		if dk == -1
			BM_table{k_+2+1,2} = BM_table{k_+2,2};
        else 		
            % 计算差值不为 -1 时，求出校正项
			correct_item = get_correct_item(dk, BM_table, field_table, m, k_);

			BM_table{k_+2+1,2} = zeros(1,length(correct_item));
            % 使用前一行的多项式进行扩展，并保持阶数一致
			tmp = [BM_table{k_+2,2},-1*ones(1,length(correct_item)-length(BM_table{k_+2,2}))];	%bao_chi wei_du yi_zhi
            
            % 更新当前多项式的每一项，利用校正项进行累加或扩展
			for i_ = 1:length(correct_item)
				if correct_item(i_) == -1 && tmp(i_) == -1
					BM_table{k_+2+1,2}(i_) = -1; % 如果两者都不存在，则保持为 -1
				elseif correct_item(i_) == -1 && tmp(i_) ~= -1
					BM_table{k_+2+1,2}(i_) = BM_table{k_+2,2}(i_); % 只存在 tmp 中的值		
				elseif correct_item(i_) ~= -1 && tmp(i_) == -1
					BM_table{k_+2+1,2}(i_) = correct_item(i_); % 只存在校正项中的值
                else 
                    % 都存在时，在有限域中执行加法操作
					BM_table{k_+2+1,2}(i_) = gf_add(correct_item(i_),tmp(i_),field_table);
				end							
			end
		end
        
        % 更新表格的其他列信息
        BM_table{k_+2+1, 1} = k_ + 1;                            % 当前迭代步
        BM_table{k_+2, 3} = dk;                                  % 差值 d(k)
        BM_table{k_+2+1, 4} = length(BM_table{k_+2+1, 2}) - 1;   % 错误定位多项式的当前阶数
        BM_table{k_+2+1, 5} = 2 * (k_ + 1) - BM_table{k_+2+1, 4}; % 2*(k+1) - 错误定位多项式的阶数

    end

    % 输出最后一行中的错误定位多项式 sigma_X
	sigma_X = BM_table{t+2,2};
end