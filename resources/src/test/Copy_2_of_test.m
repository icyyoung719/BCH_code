% 全方位测试脚本：只测试在 BCH 合理范围内的情况
clc;
fprintf('Testing main_function for k in powers of 2 and t in range 1 to 10\n');

% 打开文件写入异常信息
logFile = fopen('error_log.txt', 'w');
timeFile = fopen('time.txt', 'w');
if logFile == -1 || timeFile == -1
    error('无法打开日志文件！');
end

Time = zeros(4, 10);  % 预分配时间矩阵 (对应k_values和t的组合)

% 设置k_values为指数值（对应的k为2^k_values）
k_values = [5, 6, 7, 8];  % 指数值

% 遍历k_values
for k_idx = 1:length(k_values)
    k = 2^k_values(k_idx);  % 计算真正的k值
    
    % 遍历t的范围
    for t = 1:8
        % 计算对应的m值，满足 n = 2^m - 1 且 k = n - mt
        % 这里需要通过迭代来寻找满足条件的最小m
        m_found = false;
        for m = k_values(k_idx):12  % m的范围从2到12
            n = 2^m - 1;  % 码字长度
            if n - m * t >= k
                m_found = true;
                break;
            end
        end

        if m_found
            % 打印当前测试组合
            fprintf('Testing with k=%d, t=%d, m=%d ... ', k, t, m);
            
        %try
            % 执行主函数
            tic;
            main_function(m, t, k, logFile);
            costTime = toc;
            
            % 记录时间
            % Time(k_idx, t) = costTime;
            % fprintf('Passed (%.4f seconds)\n', costTime);
        %catch ME
            % 捕获并输出异常信息到控制台
            %fprintf('Failed\n');
            %fprintf('Error for k=%d, t=%d, m=%d: %s\n', k, t, m, ME.message);
            
            % 将错误信息写入日志文件
            %fprintf(logFile, 'Error for k=%d, t=%d, m=%d: %s\n', k, t, m, ME.message);
        %end
        else
            fprintf('No valid m found for k=%d, t=%d\n', k, t);
        end
    end

    % 输出时间矩阵的当前行到文件
    %fprintf(timeFile, '%f\t', Time(k_idx, :)); % 使用浮点数格式和制表符分隔
    %fprintf(timeFile, '\n'); % 换行
end

% 关闭文件
fclose(logFile);
fclose(timeFile);
fprintf('All tests completed.\n');
