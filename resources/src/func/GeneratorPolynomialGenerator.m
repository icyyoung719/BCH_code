function [genPoly,minpol_list,minpol_list_all] = GeneratorPolynomialGenerator(GF, m, t, primPoly)
    % 用于计算生成多项式，同时输出两种最小项多项式
    % GF: 已知的GF(2^m)有限域，由GF = de2bi(M, m, 'left-msb')产生
    % m: 域扩展的阶数，即 GF(2^m)
    % t: 纠错能力
    % primPoly: 采用的本原多项式，存放不同m对应的本原多项式的10进制表示 

    t2 = 2*t;

    % 初始化minpol_list_all
    minpol_list_all = cell(1, 2 * t); % 使用cell以支持不同长度的多项式
    p = 1; % minpol_list 的索引
    minpol_list = []; % 初始化 minpol_list
    
    % 确定域的共轭类
    coset = cosets(m, primPoly, 'nodisplay');
    
    % 生成最小多项式并填充minpol_list_all
    for idx1 = 2 : numel(coset)
        % 如果该共轭类中包含 alpha 的幂次小于 2t 的元素
        if any(find(log(coset{idx1}) < 2 * t))
            % 计算最小多项式
            tempPoly = 1;
            thisCoset = coset{idx1};
            for idx2 = 1 : length(thisCoset)
                tempPoly = conv(tempPoly, [1 thisCoset(idx2)]); % 多项式相乘
            end
            minPol = gf([zeros(1, m + 1 - length(tempPoly)), tempPoly.x], 1);
    
            % 将生成的最小多项式添加到minpol_list
            minpol_list = [minpol_list; minPol];

            % 找到minpol_list_all中第一个空位置，填充最小多项式
            i = find(cellfun(@isempty, minpol_list_all), 1); % 获取第一个空位置的索引
            if isempty(i)
                error('minpol_list_all 没有空位置，无法填充。');
            end
            % 填充到minpol_list_all
            while i <= 2 * t
                if isempty(minpol_list_all{i})
                    minpol_list_all{i} = minPol; % 复制最小多项式到minpol_list_all
                end
                i = i * 2; % 处理共轭关系
            end
    
            % 更新指针
            p = p + 1;
    
            % 如果已填充的最小多项式数超过 t，则停止
            if p > t
                break;
            end
        end
    end
    
    % 检查是否所有必要的 minpol_list_all 均已填充
    for i = 1:2*t
        if isempty(minpol_list_all{i})
            error('minpol_list_all 在索引 %d 未被正确填充，可能存在逻辑问题。', i);
        end
    end


    % 确定生成多项式
    len = size(minpol_list, 1);
    genPoly = 1;
    for i = 1:len
        genPoly = conv(genPoly, minpol_list(i,:));
    end
    % 将生成多项式的系数转换为二进制表示，并去掉前导零
    genpoly_bin = dec2bin(genPoly.x);

    % 将结果拼接起来
    decimal_number = bin2dec(genpoly_bin);
    genPoly = de2bi(decimal_number);
    
    % 去除genPoly前面的0
    firstNonZeroIdx = find(genPoly ~= 0, 1, 'first');
    if ~isempty(firstNonZeroIdx)
        genPoly = genPoly(firstNonZeroIdx:end);
    end

    % 转置
    genPoly = transpose(genPoly);
end