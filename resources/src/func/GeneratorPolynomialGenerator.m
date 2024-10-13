function [genPoly,minpol_list,minpol_list_all] = GeneratorPolynomialGenerator(GF, m, t, primPoly)
    % GF: 已知的GF(2^m)有限域，由GF = de2bi(M, m, 'left-msb')产生
    % m: 域扩展的阶数，即 GF(2^m)
    % t: 纠错能力
    % primPoly: 采用的本原多项式，存放不同m对应的本原多项式的10进制表示 

    t2 = 2*t;

    % 确定域的共轭类
    coset = cosets(m, primPoly, 'nodisplay');

    % 确定最小多项式
    minpol_list = [];
    for idx1 = 2 : numel(coset)
        if(any(find(log(coset{idx1}) < t2)))  % 如果该共轭类中包含 alpha 的幂次小于 2t 的元素
            tempPoly = 1;
            thisCoset = coset{idx1};
            for idx2 = 1 : length(thisCoset)
                tempPoly = conv(tempPoly, [1 thisCoset(idx2)]); % 多项式相乘，计算其最小项多项式
            end
            minPol = gf([zeros(1,m+1-length(tempPoly))  tempPoly.x],1);
            minpol_list = [minpol_list;minPol];
        end
    end

    % 初始化minpol_list_all
    minpol_list_all = cell(1, 2 * t); % 使用cell以支持不同长度的多项式
    filled = false(1, 2 * t); % 跟踪已填写的索引
    p=1;
    for i = 1:2*t
        % 检查 minpol_list_all[i] 是否已被填充
        if isempty(minpol_list_all{i})
            % 填充当前索引 i 及其对应的共轭类
            j = i;  % 设置当前索引为 j
            while j <= 2*t
                % 检查是否在范围内并填充
                minpol_list_all{j} = minpol_list(p, :);  % 复制最小多项式
                j = j * 2;  % 使用共轭关系（指数法则）
            end
            
            % 填充完所有共轭类后，p自增
            p = p + 1;  
            
            % 如果p > t，则结束循环
            if p > t
                break;
            end
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