% 示例：生成 GF(2^3) 的 BCH 码生成多项式
GF = gf(0:7, 3);  % 构建 GF(2^3)
m = 3;  % 域阶
t = 2;  % 纠错能力

genPoly = GeneratorPolynomialGenerator(GF, m, t);
disp('生成多项式：');
disp(genPoly);