function genPoly = GeneratorPolynomialGenerator(GF, m, t, primPoly)
    % GF: 已知的GF(2^m)有限域，由GF = de2bi(M, m, 'left-msb')产生
    % m: 域扩展的阶数，即 GF(2^m)
    % t: 纠错能力
    % primPoly: 采用的本原多项式，存放不同m对应的本原多项式的10进制表示 
    
    % 本原元素 alpha 的索引（假设 alpha 是第一个非零元素，这里用 2 作为示例）
    alpha = GF(2);  % 对应 alpha = 2
    
    % 初始化生成多项式为 1
    genPoly = gf([1], 1, primPoly);
    visitedConjugates = {};
    
    % 遍历从 alpha^0 到 alpha^(2t-1)，计算所有最小项多项式
    for i = 1:(2 * t - 1)
        % 获取 alpha^i 的共轭类
        conjClass = GetConjugateClass(alpha, i, m);
        
        % 检查是否已经遍历过这个共轭类
        if ~ismemberCell(conjClass, visitedConjugates)
            % 将共轭类添加到已遍历的列表中
            visitedConjugates = [visitedConjugates, {conjClass}];
            
            % 计算 alpha^i 的最小项多项式 phi_i(x)
            minPoly = MinimalPolynomial(conjClass, GF, m, primPoly);
            
            % 使用多项式乘法累乘最小项多项式
            genPoly = conv(genPoly, minPoly);
        end
    end
    
    % 最终生成多项式 genPoly
end

% 计算给定元素的共轭类
function conjClass = GetConjugateClass(alpha, i, m)
    conjClass = gf(zeros(1, m), m);
    for j = 0:m-1
        conjClass(j+1) = mod(i * (2^j), (2^m - 1));  % 加入每个共轭元素 alpha^(i*2^j)
    end
end

% 计算给定共轭类的最小项多项式
function minPoly = MinimalPolynomial(conjClass, GF, m, primPoly)
    % 初始化最小多项式为1
    minPoly = gf([1], m, primPoly);
    
    % 将共轭类转换成gf类的对象
    conjClassGF = gf(conjClass, m, primPoly);
    
    % 手动去重共轭类
    uniqueConjClassGF = [];
    for i = 1:length(conjClassGF)
        isUnique = true;
        for j = 1:length(uniqueConjClassGF)
            if conjClassGF(i) == uniqueConjClassGF(j)
                isUnique = false;
                break;
            end
        end
        if isUnique
            uniqueConjClassGF = [uniqueConjClassGF, conjClassGF(i)];
        end
    end
    
    % 遍历共轭类中的每个元素
    for i = 1:length(uniqueConjClassGF)
        % 使用当前共轭类元素生成的线性因子 (x - alpha^i)
        tempPoly = gf([1, -uniqueConjClassGF(i)], m, primPoly);
        
        % 多项式乘法累积（使用 gfconv 而不是 conv）
        minPoly = conv(minPoly, tempPoly);
    end
end

% 检查共轭类是否已经访问过
function result = ismemberCell(conjClass, visitedConjugates)
    result = false;
    for i = 1:length(visitedConjugates)
        if isequal(conjClass, visitedConjugates{i})
            result = true;
            return;
        end
    end
end
