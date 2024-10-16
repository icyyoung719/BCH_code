global m;
global t;
global k;

m = 5;          % 取人伽罗华域的大小GF(2^m)
t = 3;          % 纠错数
k = 16;         % 信息位数
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

% ------------------------------ Encode --------------------------------- %
% 将输入的信息转化为k位二进制串，进行编码
% 需要确保信息位长度 <= k，少于的部分会自动用0补全
input_info = 65;
if input_info >= 2^k
    error("输入的信息位过长！");
end
info = de2bi(input_info, k, 'left-msb');

% 将信息位前移，为校验位提供位置
dividend = info;
dividend(end+1:end+length(g_x)-1) = 0;

% 使用多项式除法计算校验位
checkbits = polynomial_mod(dividend, g_x);

% 信息位+校验位生成编码序列
tx_codeword = [info checkbits];
disp("编码序列：");
disp(tx_codeword);

% ------------------------------ Decode --------------------------------- %
% 模拟传输过程中的比特错误
rx_codeword = tx_codeword;
rx_codeword(4) = 1;
rx_codeword(9) = 1;
rx_codeword(22) = 0;

% 1.计算伴随式
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
    S_a(i) = poly2power(GF, polynomial_mod(generateSi_a(S_i, i), primPolyBin));
end

% 显示伴随式分量 S(x)
disp("伴随式分量S(x):");
disp(S_a);

% 2.使用berlekamp_massey算法进行译码，得到错误位置多项式

% S_a 是由之前步骤计算得到的伴随式向量
sigma_x = berlekamp_massey(S_a, t, GF, primPolyBin);

% sigma_x 就是错误定位多项式，可以用来查找根（即错误位置）
% 无法使用标准bchdec函数，因为需要编码序列位GF(2)而非GF(2^n)
