function main_function(m, t, k, logFile)
    % auto.m - 封装原main代码，接受m, t, k参数作为输入
    % 通过参数初始化变量
    
    % 设置全局变量并计算码字长度等
    n = 2^m - 1;
    p = n - k;
    d = 2 * t + 1;

    % 检验输入数据的有效性
    if k > n
        error('信息位 k 不能大于码字长度 n');
    elseif p < m * t
        error('校验位数不足以纠正 t 个错误');
    end

    % 其余代码保持不变，只需在原代码中去除m, t, k的赋值
    % 以下是main中原有的代码
    % ------------------------------ Pre-Compute --------------------------------- %
    % 由m求有限域GF(2^m)
    % primPolys存放不同m对应的本原多项式的10进制表示
    primPolys = [
        7, 11, 19, 37, 67, 131, 285, 529, 1033, 2053, 4179
    ];

    % 选择本原多项式
    primPoly = primPolys(m - 1);
    primPolyBin = de2bi(primPoly, m + 1, 'left-msb');
    M = PrimitivePolynomialGenerator(m, primPoly);
    GF = de2bi(M, m, 'left-msb');

    % 求生成多项式及最小项多项式，需要先求出各个共轭类，再求每个共轭类的最小项多项式
    [g_x , minpol_list , minpol_list_all] = GeneratorPolynomialGenerator(GF,m,t,primPoly);
    % g_x = [zeros(1, k - length(g_x)), g_x];
    % 生成field_table，是对GF的包装，用于BM译码
    field_table = cell(2^m,2);
    for i = 1:2^m
        field_table{i,1} = i-2;
        field_table{i,2} = gf(GF(i,:),1);
    end
    
    % ------------------------------ Encode --------------------------------- %
    % 将输入的信息转化为k位二进制串，进行编码
    % 需要确保信息位长度 <= k，少于的部分会自动用0补全
    
    % 随机生成一个长度小于等于 k 的二进制串
    info_length = randi(k);  % 随机生成一个 1 到 k 之间的长度
    original_info = randi([0, 1], 1, info_length);  % 生成长度为 info_length 的 0, 1 串
    % 补充 0，使得 info 的长度为 k
    info = [zeros(1, k - length(original_info)), original_info];
    
    % 将信息位前移，为校验位提供位置
    % 需要确保dividend长度为p，不足的用0补齐
    dividend = info;
    dividend(end+1:end+p) = 0;
    
    % 使用多项式除法计算校验位
    checkbits = polynomial_mod(dividend, g_x);
    checkbits = [zeros(1, p - length(checkbits)), checkbits];
    
    % 信息位+校验位生成编码序列
    tx_codeword = [info checkbits];
    
    % ------------------------- Introduce Errors----------------------------- %
    % 模拟传输过程中的比特错误
    rx_codeword = tx_codeword;
    
    err_codeword = zeros(1, n);
    % 随机选择 t 个位置，将它们设置为 1
    indices = randperm(n, t); % 随机选择 t 个不重复的位置
    err_codeword(indices) = 1;
    
    rx_codeword = bitxor(rx_codeword,err_codeword);
    

    % ------------------------------ Decode --------------------------------- %
    % Step 1.计算伴随式
    % minpol_list内存放t个最小项多项式
    % minpol_list_all 内存放2t个包含重复的最小项多项式
    
    % 初始化伴随式向量S(x)
    S_a = zeros(1, 2*t);
    
    % 对于每个 alpha^i 的伴随式，使用 alpha^i 的最小项多项式
    for i = 1:2*t
        minpol = minpol_list_all{i};
        minpol = minpol.x;
        % 计算rx_codeword与这个最小项多项式的取模，即S的2t个分量
        S_i = polynomial_mod(rx_codeword, minpol);
        % 转换为有限域中的幂次形式
        S_a(i) = poly2power(GF, polynomial_mod(generateSi_a(S_i, i, m), primPolyBin));
    end
    
    % S_a全为-1代表没有错误产生
    if(all(S_a == -1))
        disp("无错误产生");
        return;
    end
    
    % Step 2.使用berlekamp_massey算法进行译码，得到错误位置多项式
    % sigma_X 是错误定位多项式，可以用来查找根（即错误位置）
    sigma_X = berlekamp_massey(S_a,t,m,field_table);
    
    % Step 3.计算错误模式，进行纠错
    % root_sigma_X 是错误位置的根
    root_sigma_X = find_root_sigma_X(sigma_X,m,field_table);
    
    % 求出对应的错误模式err_pattern
    err_pattern = zeros(1,n);
    err_locate = mod(2^m - 1 - root_sigma_X,2^m-1);
    err_pattern(err_locate+1) = 1;
    
    err_pattern = fliplr(err_pattern);
    fx_codeword = bitxor(rx_codeword,err_pattern);
    
    % ------------------------------ Output --------------------------------- %
    % 检查译码后的结果与信息是否相同
    % disp(tx_codeword == v);
    isAllOnes = all((tx_codeword == fx_codeword) == 1);
    if isAllOnes
        disp('相同');
    else
        fprintf(logFile,'tx_codeword 与 译码结果不同');
        fprintf(logFile, 'Testing with m=%d, t=%d, k=%d ... \n', m, t, k)
        fprintf(logFile, 'tx_codeword: ');
        fprintf(logFile, '%d ', tx_codeword); % 按照格式 %d 输出每个元素
        fprintf(logFile, '\n');
        fprintf(logFile, 'er_codeword: ');
        fprintf(logFile, '%d ', err_codeword);
        fprintf(logFile, '\n');
        fprintf(logFile, 'rx_codeword: ');
        fprintf(logFile, '%d ', rx_codeword);
        fprintf(logFile, '\n');
        fprintf(logFile, 'lx_codeword: ');
        fprintf(logFile, '%d ', err_pattern);
        fprintf(logFile, '\n');
        fprintf(logFile, 'fx_codeword: ');
        fprintf(logFile, '%d ', fx_codeword);
        fprintf(logFile, '\n\n');
    end
end
