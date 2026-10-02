% =========================================================
% MATLAB 40+ 粒子轨迹同屏对比 + 速度计算 + 导出 Excel 脚本
% 参数: Pix to Micron = 0.48, FPS = 15
% 功能:
%   1. 读取 CSV 轨迹，生成真实坐标分布图与网格排列对比图
%   2. 计算各粒子全程速度指标（平均/最大/最小速度、总路程、净位移、定向性）
%   3. 计算各粒子每 5 秒（75 帧）窗口的分段平均速度
%   4. 自动生成多 Sheet Excel 报表 ('粒子速度与运动分析汇总.xlsx')
% =========================================================
clear; clc; close all;

% 1. 参数定义
pix2um = 0.48;       % 像素转微米系数 (um / pixel)
fps = 15;            % 帧率 (fps)
dt = 1 / fps;        % 帧时间间隔 (s)
interval_sec = 5;    % 分段统计时间窗口 (s)
gridSpacing = 150;   % 网格排列起始点间距 (um)
folderPath = './';   % CSV 文件所在文件夹目录

% 2. 读取并按数值顺序排序所有 CSV 文件 (1.csv, 2.csv ... 40.csv)
filePattern = fullfile(folderPath, '*.csv');
csvFiles = dir(filePattern);
if isempty(csvFiles)
    error('未在指定路径下找到 CSV 文件，请确认 folderPath 设置！');
end

[~, sortIdx] = sort(cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), {csvFiles.name}));
csvFiles = csvFiles(sortIdx);
numFiles = length(csvFiles);
fprintf('成功找到 %d 个轨迹文件，开始读取数据与计算速度...\n', numFiles);

numCols = ceil(sqrt(numFiles));
numRows = ceil(numFiles / numCols);

% 3. 初始化绘图窗口
figure('Name', '多粒子轨迹网格排列与速度分析汇总', 'Color', 'w', 'Position', [100, 100, 1200, 520]);
colors = turbo(numFiles);

% --- 子图 1：原始物理坐标轨迹 ---
subplot(1, 2, 1); hold on; grid on; box on;
title('所有粒子的真实物理位置分布 (原图)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('X (\mum)', 'FontSize', 11); ylabel('Y (\mum)', 'FontSize', 11); axis equal;

% --- 子图 2：网格平移展示 ---
subplot(1, 2, 2); hold on; grid on; box on;
title(['粒子的网格排列轨迹（间距 ' num2str(gridSpacing) ' \mum）'], 'FontSize', 12, 'FontWeight', 'bold');
xlabel('X (\mum)', 'FontSize', 11); ylabel('Y (\mum)', 'FontSize', 11); axis equal;

% 4. 数据统计容器初始化
summaryData = cell(numFiles, 9);         % 总体统计表
framesPerInterval = interval_sec * fps; % 5s 对应的帧步数 (75 步)
particleSpeeds = cell(numFiles, 1);     % 保存每个粒子的逐帧速度
maxIntervals = 0;
maxFrames = 0;

% 5. 循环读取 CSV 文件、绘制图形并计算速度
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

    % --- 轨迹绘制 (子图 1: 真实位置) ---
    subplot(1, 2, 1);
    plot(x_um, y_um, '-', 'LineWidth', 1.2, 'Color', colors(k, :), 'HandleVisibility', 'off');
    plot(x_um(1), y_um(1), 'o', 'MarkerSize', 3.5, 'MarkerFaceColor', colors(k, :), 'MarkerEdgeColor', 'none');
    plot(x_um(end), y_um(end), 's', 'MarkerSize', 4.5, 'MarkerFaceColor', colors(k, :), 'MarkerEdgeColor', 'k');

    % --- 轨迹绘制 (子图 2: 网格平移) ---
    subplot(1, 2, 2);
    rowIdx = floor((k-1) / numCols); colIdx = mod(k-1, numCols);
    x_offset = colIdx * gridSpacing; y_offset = -rowIdx * gridSpacing;
    x_grid = (x_um - x_um(1)) + x_offset;
    y_grid = (y_um - y_um(1)) + y_offset;

    plot(x_grid, y_grid, '-', 'LineWidth', 1.2, 'Color', colors(k, :), 'HandleVisibility', 'off');
    plot(x_grid(1), y_grid(1), 'o', 'MarkerSize', 3.5, 'MarkerFaceColor', colors(k, :), 'MarkerEdgeColor', 'none');
    plot(x_grid(end), y_grid(end), 's', 'MarkerSize', 4.5, 'MarkerFaceColor', colors(k, :), 'MarkerEdgeColor', 'k');
    if numFiles <= 60
        text(x_grid(1), y_grid(1) + gridSpacing*0.1, ['p=' num2str(k)], 'Color', 'k', ...
            'FontSize', 8, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end

    % -------------------------------------------------------------
    % 速度与运动参数计算
    % -------------------------------------------------------------
    dx = diff(x_um);
    dy = diff(y_um);
    step_dist = sqrt(dx.^2 + dy.^2);  % 逐帧位移 (um)
    step_speed = step_dist / dt;      % 逐帧瞬时速度 (um/s)
    particleSpeeds{k} = step_speed;

    total_time = (N - 1) * dt;                  % 总时间 (s)
    total_dist = sum(step_dist);                % 累积总路程 (um)
    net_displ = sqrt((x_um(end)-x_um(1))^2 + (y_um(end)-y_um(1))^2); % 净位移 (um)
    mean_speed = mean(step_speed);              % 平均速度 (um/s)
    max_speed = max(step_speed);                % 最大速度 (um/s)
    min_speed = min(step_speed);                % 最小速度 (um/s)
    directionality = net_displ / total_dist;   % 定向性指数 (D/S)

    % 保存总体数据
    particleID = sprintf('Particle_%02d', k);
    summaryData(k, :) = {particleID, N, total_time, total_dist, net_displ, mean_speed, max_speed, min_speed, directionality};

    % 记录最多的 5s 间隔数量
    numIntervals = floor(length(step_speed) / framesPerInterval);
    if numIntervals > maxIntervals
        maxIntervals = numIntervals;
    end
end

% 图像格式修饰
colormap(turbo(numFiles));
c = colorbar; c.Label.String = '粒子编号 (文件索引)'; caxis([1 numFiles]);

% -------------------------------------------------------------
% 6. 构建每 5s 分段平均速度表
% -------------------------------------------------------------
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
            intervalData{k, i+1} = mean(v(idx_start:idx_end));
        else
            intervalData{k, i+1} = NaN;
        end
    end
    intervalData{k, end} = mean(v);
end

% -------------------------------------------------------------
% 7. 构建逐帧实时速度明细表
% -------------------------------------------------------------
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

% -------------------------------------------------------------
% 8. 自动导出数据至 Excel 文件 ('粒子速度与运动分析汇总.xlsx')
% -------------------------------------------------------------
excelFileName = '粒子速度与运动分析汇总.xlsx';

summaryHeaders = {'粒子编号', '总帧数', '总时长_s', '总路程_um', '净位移_um', '平均速度_um_s', '最大速度_um_s', '最小速度_um_s', '定向性指数'};

T_summary = cell2table(summaryData, 'VariableNames', summaryHeaders);
T_interval = cell2table(intervalData, 'VariableNames', intervalHeaders);
T_detail = cell2table(detailData, 'VariableNames', detailHeaders);

% 写入同一个 Excel 的三个 Sheet 页
writetable(T_summary, excelFileName, 'Sheet', '总体速度统计');
writetable(T_interval, excelFileName, 'Sheet', '每5秒分段平均速度');
writetable(T_detail, excelFileName, 'Sheet', '逐帧实时速度明细');

fprintf('\n处理完毕！\n1. 轨迹网格对比图已生成；\n2. 速度统计结果已保存至：%s\n', fullfile(pwd, excelFileName));