% 全方位测试脚本：只测试在 BCH 合理范围内的情况
fprintf('Testing main_function for m in range 2 to 12\n');
max_m = 12;

% 打开文件写入异常信息
logFile = fopen('error_log.txt', 'w');
if logFile == -1
    error('无法打开日志文件！');
end

% 根据 m 的值调整最大 t 的范围
for m = 5:max_m
    n = 2^m - 1;  % 码字长度
    max_t = floor((n - m) / (2 * m));  % 根据 m 计算最大 t 的值

    % t 的合理范围
    for t = 1:max_t
        % 对于每个 t，确定 k 的合理范围
        k_min = m * t;  % 最小信息位数 k

        for k = k_min:(n - m * t)  % k 的范围必须满足最低限制
            % 打印当前测试组合
            fprintf('Testing with m=%d, t=%d, k=%d ... ', m, t, k);
            
            try
                % 执行主函数
                main_function(m, t, k, logFile);
                
                % 确认无错误产生
                % fprintf('Passed\n');
            catch ME
                % 捕获并输出异常信息到控制台
                fprintf('Failed\n');
                fprintf('Error for m=%d, t=%d, k=%d: %s\n', m, t, k, ME.message);
                
                % 将错误信息写入日志文件
                fprintf(logFile, 'Error for m=%d, t=%d, k=%d: %s\n', m, t, k, ME.message);
            end
        end
    end
end

% 关闭文件
fclose(logFile);

fprintf('All tests completed.\n');
