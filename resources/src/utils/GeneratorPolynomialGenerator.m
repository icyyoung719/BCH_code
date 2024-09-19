function [genpoly, t] = GeneratorPolynomialGenerator(GF, m, t, primPoly)
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
        if(any(find(log(coset{idx1}) < t2)))  % coset contains a power of alpha < 2t 
            tempPoly = 1;
            thisCoset = coset{idx1};
            for idx2 = 1 : length(thisCoset)
                tempPoly = conv(tempPoly, [1 thisCoset(idx2)]);
            end
            minPol = gf([zeros(1,m+1-length(tempPoly))  tempPoly.x],1);
            minpol_list = [minpol_list;minPol];
        end
    end

    % 确定生成多项式
    len = size(minpol_list, 1);
    genpoly = 1;
    for i = 1:len
        genpoly = conv(genpoly, minpol_list(i,:));
    end
    genpoly = genpoly(end-(GF-m) :end);
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
