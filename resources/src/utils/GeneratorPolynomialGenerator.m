function genPoly = GeneratorPolynomialGenerator(GF, m, t, primPoly)
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
