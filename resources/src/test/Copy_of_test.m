% 全方位测试脚本：只测试在 BCH 合理范围内的情况
% 测试m：1到11，选择最大的k
fprintf('Testing main_function for m in range 2 to 12\n');
max_m = 12;

% 打开文件写入异常信息
logFile = fopen('error_log.txt', 'w');
timeFile = fopen('time.txt','w');
if logFile == -1
    error('无法打开日志文件！');
end

Time = zeros(max_m,max_m);

% 根据 m 的值调整最大 t 的范围
for m = 1:11
    n = 2^m - 1;  % 码字长度
    max_t = floor((n - m) / (2 * m));  % 根据 m 计算最大 t 的值
    if max_t > m
        max_t = m;
    end

    % t 的合理范围
    for t = 1:max_t
        % 对于每个 t，确定 k 的合理范围
        k = n - m*t;  % 最小信息位数 k

        if k > 0 
            % 打印当前测试组合
            fprintf('Testing with m=%d, t=%d, k=%d ... ', m, t, k);
            
            try
                % 执行主函数
                tic;
                main_function(m, t, k, logFile);
                costTime = toc;
                Time(m,t) = costTime;
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
    fprintf(timeFile, '%f\t', Time(m, :)); % 使用浮点数格式和制表符分隔
    fprintf(timeFile, '\n'); % 换行
end

% 关闭文件
fclose(logFile);
%writematrix(Time,timeFile);
%for i = 1:size(Time, 1)

%end
%disp(Time);
fclose(timeFile);
fprintf('All tests completed.\n');
