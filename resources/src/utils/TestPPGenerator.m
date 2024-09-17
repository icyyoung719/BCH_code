% 从命令行读取用户输入
fprintf('请输入 m 的值 (整数): ');
m = input('');

% 生成 M 数组
M = PrimitivePolynomialGenerator(m);

% 显示结果
disp('生成的 M 数组:');
disp(M);