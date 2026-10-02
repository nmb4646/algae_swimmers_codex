% =========================================================
% MATLAB 粒子速度可视化 + 每 3 秒分段速度分析 + Excel 导出脚本
% 参数: Pix to Micron = 0.48, FPS = 15, 时间窗口 = 3 秒
% 功能:
%   1. 读取 CSV 轨迹，计算逐帧瞬时速度与每 3 秒分段平均速度
%   2. 界面画图展示：直接看粒子速度对比及【速度随时间下降（变慢）的趋势折线】
%   3. 导出 Excel 数据表（包含 3 秒分段速度、总体速度、逐帧明细）
% =========================================================
clear; clc; close all;

% 1. 参数定义
pix2um = 0.48;       % 像素转微米系数 (um / pixel)
fps = 15;            % 帧率 (fps)
dt = 1 / fps;        % 帧时间间隔 (s)
interval_sec = 3;    % 分段统计时间窗口 (3 秒)
folderPath = './';   % CSV 文件所在文件夹

% 2. 获取并按阿拉伯数字顺序排序所有 CSV 文件 (1.csv, 2.csv ... 40.csv)
filePattern = fullfile(folderPath, '*.csv');
csvFiles = dir(filePattern);
if isempty(csvFiles)
    error('未在指定路径下找到 CSV 文件，请确认 folderPath 设置！');
end

[~, sortIdx] = sort(cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), {csvFiles.name}));
csvFiles = csvFiles(sortIdx);
numFiles = length(csvFiles);
fprintf('成功找到 %d 个轨迹文件，正在分析速度...\n', numFiles);

% 3. 数据容器初始化
summaryData = cell(numFiles, 9);         % 总体统计表
framesPerInterval = interval_sec * fps; % 3s 对应的帧步数 (15 * 3 = 45 步)
particleSpeeds = cell(numFiles, 1);     % 保存每个粒子的逐帧速度
maxIntervals = 0;
maxFrames = 0;

% 4. 循环读取 CSV 文件并计算速度
for k = 1:numFiles
    filePath = fullfile(folderPath, csvFiles(k).name);
    opts = detectImportOptions(filePath);
    opts.VariableNamingRule = 'preserve';
    tbl = readtable(filePath, opts);

    % 自动提取 X 和 Y 列
    colNames = tbl.Properties.VariableNames;
    xIdx = find(strcmpi(strtrim(colNames), 'X'), 1);
    yIdx = find(strcmpi(strtrim(colNames), 'Y'), 1);

    if isempty(xIdx) || isempty(yIdx)
        x_pix = tbl{:, 3}; y_pix = tbl{:, 4};
    else
        x_pix = tbl{:, xIdx}; y_pix = tbl{:, yIdx};
    end

    x_um = x_pix * pix2um;
    y_um = y_pix * pix2um;
    N = length(x_um);
    if N > maxFrames, maxFrames = N; end

    % 计算逐帧位移与瞬时速度
    dx = diff(x_um);
    dy = diff(y_um);
    step_dist = sqrt(dx.^2 + dy.^2);  % 逐帧位移 (um)
    step_speed = step_dist / dt;      % 逐帧速度 (um/s)
    particleSpeeds{k} = step_speed;

    total_time = (N - 1) * dt;
    total_dist = sum(step_dist);
    net_displ = sqrt((x_um(end)-x_um(1))^2 + (y_um(end)-y_um(1))^2);
    mean_speed = mean(step_speed);
    max_speed = max(step_speed);
    min_speed = min(step_speed);
    directionality = net_displ / total_dist;

    particleID = sprintf('Particle_%02d', k);
    summaryData(k, :) = {particleID, N, total_time, total_dist, net_displ, mean_speed, max_speed, min_speed, directionality};

    numIntervals = floor(length(step_speed) / framesPerInterval);
    if numIntervals > maxIntervals
        maxIntervals = numIntervals;
    end
end

% 5. 构建每 3s 分段平均速度矩阵（用于画图与导出 Excel）
intervalMatrix = NaN(numFiles, maxIntervals); % 存储各时间段数值用于绘图
intervalData = cell(numFiles, maxIntervals + 2);
intervalHeaders = cell(1, maxIntervals + 2);
intervalHeaders{1} = '粒子编号';

for i = 1:maxIntervals
    intervalHeaders{i+1} = sprintf('%ds_%ds_平均速度', (i-1)*interval_sec, i*interval_sec);
end
intervalHeaders{end} = '全程平均速度';

for k = 1:numFiles
    v = particleSpeeds{k};
    intervalData{k, 1} = sprintf('Particle_%02d', k);

    for i = 1:maxIntervals
        idx_start = (i - 1) * framesPerInterval + 1;
        idx_end = i * framesPerInterval;

        if idx_end <= length(v)
            seg_mean = mean(v(idx_start:idx_end));
            intervalMatrix(k, i) = seg_mean;
            intervalData{k, i+1} = seg_mean;
        end
    end
    intervalData{k, end} = mean(v);
end

% 6. 【重点修改：在 MATLAB 窗口中直接显示速度分析图表】
figure('Name', '粒子速度与趋势分析看板', 'Color', 'w', 'Position', [100, 100, 1200, 520]);
colors = turbo(numFiles);

% --- 子图 1：各粒子全程平均速度对比 (柱状图) ---
subplot(1, 2, 1);
allMeanSpeeds = [summaryData{:, 6}]; % 提取全程平均速度
b = bar(1:numFiles, allMeanSpeeds, 'FaceColor', 'flat');
b.CData = colors;
grid on; box on;
title('各粒子全程平均速度对比', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('粒子编号 (Index)', 'FontSize', 11);
ylabel('平均速度 (\mum/s)', 'FontSize', 11);
xlim([0, numFiles + 1]);

% --- 子图 2：各粒子每 3 秒分段速度变化趋势 (折线图) ---
subplot(1, 2, 2); hold on; grid on; box on;
timeX = (1:maxIntervals) * interval_sec; % X 轴时间节点 (3s, 6s, 9s...)
for k = 1:numFiles
    plot(timeX, intervalMatrix(k, :), '-o', 'LineWidth', 1.2, ...
        'Color', colors(k, :), 'MarkerSize', 4, 'DisplayName', sprintf('P%02d', k));
end
title('各粒子每 3 秒速度变化趋势（线往下滑即表示减速）', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('时间点 (s)', 'FontSize', 11);
ylabel('分段平均速度 (\mum/s)', 'FontSize', 11);
colormap(turbo(numFiles));
c = colorbar; c.Label.String = '粒子编号'; caxis([1 numFiles]);

% 7. 构建逐帧实时速度明细表
detailData = cell(maxFrames - 1, numFiles + 2);
detailHeaders = cell(1, numFiles + 2);
detailHeaders{1} = '帧号_Frame';
detailHeaders{2} = '时间_s';
for k = 1:numFiles
    detailHeaders{k+2} = sprintf('P%02d_速度_um_s', k);
end

for f = 1:(maxFrames - 1)
    detailData{f, 1} = f;
    detailData{f, 2} = f * dt;
    for k = 1:numFiles
        v = particleSpeeds{k};
        if f <= length(v)
            detailData{f, k+2} = v(f);
        else
            detailData{f, k+2} = NaN;
        end
    end
end

% 8. 导出数据至 Excel 文件 ('粒子速度与运动分析汇总.xlsx')
excelFileName = '粒子速度与运动分析汇总.xlsx';

summaryHeaders = {'粒子编号', '总帧数', '总时长_s', '总路程_um', '净位移_um', '平均速度_um_s', '最大速度_um_s', '最小速度_um_s', '定向性指数'};

T_summary = cell2table(summaryData, 'VariableNames', summaryHeaders);
T_interval = cell2table(intervalData, 'VariableNames', intervalHeaders);
T_detail = cell2table(detailData, 'VariableNames', detailHeaders);

writetable(T_summary, excelFileName, 'Sheet', '总体速度统计');
writetable(T_interval, excelFileName, 'Sheet', '每3秒分段平均速度');
writetable(T_detail, excelFileName, 'Sheet', '逐帧实时速度明细');

fprintf('\n运行完成！\n1. MATLAB 界面已弹出速度分析与减速趋势折线图；\n2. 详细数据已保存至 Excel：%s\n', fullfile(pwd, excelFileName));