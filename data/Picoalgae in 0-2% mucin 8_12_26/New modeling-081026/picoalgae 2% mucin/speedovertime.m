% =========================================================
% MATLAB 粒子速度可视化 + 每 3 秒分段速度趋势 + 文件名数字精确对应 + Excel 导出
% 参数: Pix to Micron = 0.48, FPS = 15, 时间窗口 = 3 秒
% 特性: 图像编号 (P1, P2...) 严格对应原始 CSV 文件名 (1.csv, 2.csv...)
% =========================================================
clear; clc; close all;

% 1. 参数定义
pix2um = 0.48;       % 像素转微米系数 (um / pixel)
fps = 15;            % 帧率 (fps)
dt = 1 / fps;        % 帧时间间隔 (s)
interval_sec = 3;    % 分段统计时间窗口 (3 秒)
folderPath = './';   % CSV 文件所在文件夹

% 2. 获取 CSV 文件并提取文件名中的阿拉伯数字进行精确排序
filePattern = fullfile(folderPath, '*.csv');
csvFiles = dir(filePattern);
if isempty(csvFiles)
    error('未在指定路径下找到 CSV 文件，请确认 folderPath 设置！');
end

% 提取文件名中的数字并升序排列
fileNumbers = cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), {csvFiles.name});
[fileNumbers, sortIdx] = sort(fileNumbers);
csvFiles = csvFiles(sortIdx);
numFiles = length(csvFiles);

fprintf('成功找到 %d 个轨迹文件，正在按文件名数字关联分析速度...\n', numFiles);

% 3. 数据容器初始化
summaryData = cell(numFiles, 9);         % 总体统计表
framesPerInterval = interval_sec * fps; % 3s 对应的帧步数 (15 * 3 = 45 步)
particleSpeeds = cell(numFiles, 1);     % 保存每个粒子的逐帧速度
maxIntervals = 0;
maxFrames = 0;

% 4. 循环读取 CSV 文件并计算速度
for k = 1:numFiles
    pNum = fileNumbers(k); % 当前文件对应的真实阿拉伯数字 (如 1, 2, 12 等)
    
    filePath = fullfile(folderPath, csvFiles(k).name);
    opts = detectImportOptions(filePath);
    opts.VariableNamingRule = 'preserve';
    tbl = readtable(filePath, opts);

    % --- 自动解析 X, Y 坐标列 ---
    colNames = tbl.Properties.VariableNames;
    numTableCols = width(tbl);
    
    xIdx = find(~cellfun(@isempty, regexpi(colNames, '^x$|^x_|^x\s|\(x\)|pos.*x')), 1);
    yIdx = find(~cellfun(@isempty, regexpi(colNames, '^y$|^y_|^y\s|\(y\)|pos.*y')), 1);
    
    if isempty(xIdx), xIdx = find(contains(upper(colNames), 'X'), 1); end
    if isempty(yIdx), yIdx = find(contains(upper(colNames), 'Y'), 1); end

    if isempty(xIdx) || isempty(yIdx)
        if numTableCols == 2
            xIdx = 1; yIdx = 2;
        elseif numTableCols >= 4
            xIdx = 3; yIdx = 4;
        elseif numTableCols == 3
            xIdx = 1; yIdx = 2;
        else
            error('文件 %s 列数不符合要求，无法读取 X, Y！', csvFiles(k).name);
        end
    end

    x_raw = tbl{:, xIdx};
    y_raw = tbl{:, yIdx};
    if iscell(x_raw), x_raw = str2double(x_raw); end
    if iscell(y_raw), y_raw = str2double(y_raw); end

    x_um = x_raw * pix2um;
    y_um = y_raw * pix2um;
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
    mean_speed = mean(step_speed, 'omitnan');
    max_speed = max(step_speed, [], 'omitnan');
    min_speed = min(step_speed, [], 'omitnan');
    directionality = net_displ / total_dist;

    % 使用与文件名对应的真实数字进行编号 (例: Particle_12)
    particleID = sprintf('Particle_%d', pNum);
    summaryData(k, :) = {particleID, N, total_time, total_dist, net_displ, mean_speed, max_speed, min_speed, directionality};

    numIntervals = floor(length(step_speed) / framesPerInterval);
    if numIntervals > maxIntervals
        maxIntervals = numIntervals;
    end
end

% 5. 构建每 3s 分段平均速度矩阵
intervalMatrix = NaN(numFiles, maxIntervals);
intervalData = cell(numFiles, maxIntervals + 2);
intervalHeaders = cell(1, maxIntervals + 2);
intervalHeaders{1} = '粒子编号';

for i = 1:maxIntervals
    intervalHeaders{i+1} = sprintf('%ds_%ds_平均速度', (i-1)*interval_sec, i*interval_sec);
end
intervalHeaders{end} = '全程平均速度';

for k = 1:numFiles
    pNum = fileNumbers(k);
    v = particleSpeeds{k};
    intervalData{k, 1} = sprintf('Particle_%d', pNum);

    for i = 1:maxIntervals
        idx_start = (i - 1) * framesPerInterval + 1;
        idx_end = i * framesPerInterval;

        if idx_end <= length(v)
            seg_mean = mean(v(idx_start:idx_end), 'omitnan');
            intervalMatrix(k, i) = seg_mean;
            intervalData{k, i+1} = seg_mean;
        end
    end
    intervalData{k, end} = mean(v, 'omitnan');
end

% 6. 【画图展示：精确显示与文件名对应的阿拉伯数字】
figure('Name', '粒子速度与趋势分析看板', 'Color', 'w', 'Position', [80, 80, 1300, 560]);
colors = turbo(numFiles);

% --- 子图 1：各粒子全程平均速度对比 ---
subplot(1, 2, 1);
allMeanSpeeds = [summaryData{:, 6}]; % 提取全程平均速度
b = bar(1:numFiles, allMeanSpeeds, 'FaceColor', 'flat');
b.CData = colors;
grid on; box on;
title('各粒子全程平均速度对比', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('文件名对应粒子编号', 'FontSize', 11);
ylabel('平均速度 (\mum/s)', 'FontSize', 11);
xlim([0, numFiles + 1]);

% 设置 X 轴刻度显示为精确的文件名数字 P1, P2... P44
xticks(1:numFiles);
xlabels = arrayfun(@(x) sprintf('P%d', fileNumbers(x)), 1:numFiles, 'UniformOutput', false);
xticklabels(xlabels);
xtickangle(45);
set(gca, 'FontSize', 8);

% 柱子上方标注与文件名完全相符的数字
maxBarVal = max(allMeanSpeeds);
for k = 1:numFiles
    text(k, allMeanSpeeds(k) + maxBarVal*0.02, sprintf('P%d', fileNumbers(k)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 7.5, 'Rotation', 90, 'FontWeight', 'bold', 'Color', colors(k, :));
end
ylim([0, maxBarVal * 1.25]);

% --- 子图 2：各粒子每 3 秒速度变化趋势 (线末端直接标注 P数字) ---
subplot(1, 2, 2); hold on; grid on; box on;
timeX = (1:maxIntervals) * interval_sec; % X 轴时间节点 (3s, 6s, 9s...)
for k = 1:numFiles
    pNum = fileNumbers(k);
    plot(timeX, intervalMatrix(k, :), '-o', 'LineWidth', 1.2, ...
        'Color', colors(k, :), 'MarkerSize', 4, 'DisplayName', sprintf('P%d', pNum));
    
    % 在折线末端直接标注对应文件名的数字
    valid_idx = find(~isnan(intervalMatrix(k, :)), 1, 'last');
    if ~isempty(valid_idx)
        text(timeX(valid_idx) + 0.1, intervalMatrix(k, valid_idx), sprintf(' P%d', pNum), ...
            'Color', colors(k, :), 'FontSize', 8, 'FontWeight', 'bold', ...
            'VerticalAlignment', 'middle');
    end
end

title('各粒子每 3 秒速度变化趋势（末端标注文件名数字）', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('时间点 (s)', 'FontSize', 11);
ylabel('分段平均速度 (\mum/s)', 'FontSize', 11);
xlim([min(timeX)-0.5, max(timeX) + interval_sec*0.6]);
colormap(turbo(numFiles));
c = colorbar; c.Label.String = '粒子文件名序号'; caxis([min(fileNumbers) max(fileNumbers)]);

% 7. 构建逐帧实时速度明细表 (表头同样与文件名匹配)
detailData = cell(maxFrames - 1, numFiles + 2);
detailHeaders = cell(1, numFiles + 2);
detailHeaders{1} = '帧号_Frame';
detailHeaders{2} = '时间_s';
for k = 1:numFiles
    detailHeaders{k+2} = sprintf('P%d_速度_um_s', fileNumbers(k));
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

fprintf('\n运行完成！\n1. 图表与 Excel 中的 P1, P2, P12... 已与 CSV 文件名中的阿拉伯数字 1:1 精确绑定；\n2. 详细数据已导出至：%s\n', fullfile(pwd, excelFileName));