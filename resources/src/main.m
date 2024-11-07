% Found Errors:
% m=10 , t=3，2
% m=7,t =4
% m=7 t=3 结果可能正确，也可能出错
% ------------------------------ Begin --------------------------------- %
global m;
global t;
global k;

m = 10;          % 取人伽罗华域的大小GF(2^m)
t = 2;          % 纠错数
k = 7;         % 信息位数
n = 2^m - 1;    % 编码后长度
p = n - k;      % 检验位数
d = 2*t + 1;    % 最小汉明距离
% 检验输入的数据是否冲突
if k > n
    error('信息位 k 不能大于码字长度 n');
elseif p < m*t
    error('校验位数不足以纠正 t 个错误');
end

% 由m求有限域GF(2^m)
% primPolys存放不同m对应的本原多项式的10进制表示 
primPolys = [
    7,      % 对应m=2，X^2 + X + 1 -> 111 (binary) = 7 (decimal)
    11,     % 对应m=3，X^3 + X + 1 -> 1011 (binary) = 11 (decimal)
    19,     % 对应m=4，X^4 + X + 1 -> 10011 (binary) = 19 (decimal)
    37,     % 对应m=5，X^5 + X^2 + 1 -> 100101 (binary) = 37 (decimal)
    67,     % 对应m=6，X^6 + X + 1 -> 1000011 (binary) = 67 (decimal)
    131,    % 对应m=7，X^7 + X + 1 -> 10000011 (binary) = 131 (decimal)
    285,    % 对应m=8，X^8 + X^4 + X^3 + X^2 + 1 -> 100011101 (binary) = 285 (decimal)
    529,    % 对应m=9，X^9 + X^4 + 1 -> 1000010001 (binary) = 529 (decimal)
    1033,   % 对应m=10，X^10 + X^3 + 1 -> 10000001001 (binary) = 1033 (decimal)
    2053,   % 对应m=11，X^11 + X^2 + 1 -> 100000000101 (binary) = 2053 (decimal)
    4179    % 对应m=12，X^12 + X^6 + X^4 + X + 1 -> 100000100011 (binary) = 4179 (decimal)
];

% 由m选择预先准备的本原多项式
primPoly = primPolys(m - 1);  % 根据m选择相应的本原多项式
primPolyBin = de2bi(primPoly, m + 1, 'left-msb');  % 本原多项式的二进制形式
M = PrimitivePolynomialGenerator(m,primPoly);
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

% 随机获取长度小于等于k的信息数据
max_value = 2^k - 1;
% input_info = randi(max_value);
input_info = 65;
info_length = floor(log2(input_info)) + 1;
if info_length > k
    error("输入的信息位过长，无法用 k 位二进制数表示！");
end

info = de2bi(input_info, k, 'left-msb');

% 将信息位前移，为校验位提供位置
% 需要确保dividend长度为p，不足的用0补齐
dividend = info;
dividend(end+1:end+p) = 0;

% 使用多项式除法计算校验位
checkbits = polynomial_mod(dividend, g_x);
checkbits = [zeros(1, p - length(checkbits)), checkbits];

% 信息位+校验位生成编码序列
tx_codeword = [info checkbits];
% disp("编码序列：");
% disp(tx_codeword);

% ------------------------- Introduce Errors----------------------------- %
% 模拟传输过程中的比特错误
rx_codeword = tx_codeword;
%{
err_codeword = zeros(1, n);
% 随机选择 t 个位置，将它们设置为 1
indices = randperm(n, t); % 随机选择 t 个不重复的位置
err_codeword(indices) = 1;

rx_codeword = bitxor(rx_codeword,err_codeword);
%}
rx_codeword(4) = 1;
rx_codeword(9) = 1;
rx_codeword(22) = 0;

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

% 显示伴随式分量 S(x)
disp("伴随式分量S_a:");
disp(S_a);

% S_a全为-1代表没有错误产生
if(all(S_a == -1))
    disp("无错误产生");
    return;
end

% Step 2.使用berlekamp_massey算法进行译码，得到错误位置多项式
% sigma_X 是错误定位多项式，可以用来查找根（即错误位置）
sigma_X = berlekamp_massey(S_a,t,m,field_table);
disp("错误位置多项式:");
disp(sigma_X);

% Step 3.计算错误模式，进行纠错
% root_sigma_X 是错误位置的根
root_sigma_X = find_root_sigma_X(sigma_X,m,field_table);

% 求出对应的错误模式err_pattern
err_pattern = zeros(1,n);
err_locate = mod(2^m - 1 - root_sigma_X,2^m-1);
err_pattern(err_locate+1) = 1;

err_pattern = fliplr(err_pattern);
fx_codeword = bitxor(rx_codeword,err_pattern);
% disp("错误位置：");
% disp(err_pattern);

% ------------------------------ Output --------------------------------- %
% 检查译码后的结果与信息是否相同
% disp(tx_codeword == v);
isAllOnes = all((tx_codeword == fx_codeword) == 1);
if isAllOnes
    disp('tx_codeword 与 译码结果相同');
else
    disp('tx_codeword 与 译码结果不同');
end

