function [decoded, errors] = bch_decode(r, n, k, t)
% r: 接收到的码字（含噪声）
% n: 码长
% k: 信息位长度
% t: 可以纠正的错误数

% 计算生成多项式的根
m = log2(n+1); % GF(2^m) 中的元素
primpoly = gfprimdf(m, 'min'); % 获取本原多项式
alpha = gf(2, m, primpoly); % 本原元素
g = bchgenpoly(n, k); % 生成BCH码的生成多项式

% 计算伴随式
S = zeros(1, 2*t);
for i = 1:2*t
    S(i) = polyval(gf(r, m, primpoly), alpha^i);
end

% Berlekamp-Massey算法
L = 0; % 错误定位器多项式的阶
Lambda = [1, zeros(1, 2*t)]; % 初始化错误定位器多项式
B = [1, zeros(1, 2*t)]; % 辅助多项式
Omega = zeros(1, 2*t); % 错误评估器多项式

for i = 1:2*t
    Delta = S(i);
    for j = 1:L
        Delta = Delta + Lambda(j+1) * S(i-j);
    end
    
    if Delta ~= 0
        T = Lambda;
        Lambda = Lambda - Delta * fliplr(B) .* (alpha.^(0:i-1));
        
        if 2*L <= i
            L = i - L;
            B = T / Delta;
            B = [B, zeros(1, 2*t-length(B))];
        else
            B = B * alpha^(-Delta);
        end
    end
end

% 找到错误位置
errors = findroots(Lambda, m, primpoly);

% 纠正错误
decoded = r;
for e = 1:length(errors)
    idx = mod(n - errors(e), n);
    decoded(idx+1) = decoded(idx+1) + 1; % 假设二进制编码，简单翻转位
end

% 返回译码后的数据
decoded = decoded(1:k);
end

function roots = findroots(poly, m, primpoly)
% 寻找多项式在GF(2^m)中的根
gf_poly = gf(poly, m, primpoly);
roots = [];
for i = 1:2^m-1
    if polyval(gf_poly, gf(i, m, primpoly)) == 0
        roots = [roots, i];
    end
end
end