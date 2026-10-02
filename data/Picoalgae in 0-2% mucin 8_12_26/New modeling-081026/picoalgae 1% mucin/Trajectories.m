% =========================================================
% MATLAB 40+ 粒子轨迹网格排列 + 瞬时速度伪彩映射 (固定 Colorbar 20-80 um/s)
% 参数: Pix to Micron = 0.48, FPS = 15
% 功能: 轨迹线条按瞬时速度(um/s)渐变着色，Colorbar固定范围 20~80 um/s
% =========================================================

clear; clc; close all;

% 1. 参数定义
pix2um = 0.48;       % 像素转微米系数 (um / pixel)
fps = 15;            % 帧率 (fps)
dt = 1 / fps;        % 帧时间间隔 (s)
folderPath = './';   % CSV 文件所在文件夹 ('./' 表示当前脚本同级目录)
gridSpacing = 150;   % 微米 (um), 网格排列时起始点之间的间距

% 色彩范围限制 (固定为 20 - 80 um/s)
minSpeedLimit = 20;  % 最低速度 (对应冷色/蓝)
maxSpeedLimit = 80;  % 最高速度 (对应暖色/红)

% 2. 获取并按阿拉伯数字顺序排序所有 CSV 文件 (1.csv, 2.csv ... 40.csv)
filePattern = fullfile(folderPath, '*.csv');
csvFiles = dir(filePattern);

if isempty(csvFiles)
    error('未在指定路径下找到 CSV 文件，请确认 folderPath 设置！');
end

% 提取文件名中的阿拉伯数字并升序排列
[~, sortIdx] = sort(cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), {csvFiles.name}));
csvFiles = csvFiles(sortIdx);

numFiles = length(csvFiles);
fprintf('成功找到 %d 个轨迹文件，开始计算速度并绘制...\n', numFiles);

% 计算网格行列数
numCols = ceil(sqrt(numFiles)); 
numRows = ceil(numFiles / numCols); 

% 3. 初始化窗口与画板
figure('Name', '粒子轨迹 + 瞬时速度伪彩映射 (20-80 um/s)', 'Color', 'w', 'Position', [100, 100, 1250, 550]);

% 选择速度伪彩色板 (如 turbo, jet, viridis, hot)
currentColormap = turbo;

% --- 子图 1：原始物理坐标轨迹 ---
subplot(1, 2, 1);
hold on; grid on; box on;
title('真实物理位置轨迹 (颜色表示瞬时速度)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('X (\mum)', 'FontSize', 11);
ylabel('Y (\mum)', 'FontSize', 11);
axis equal;

% --- 子图 2：网格平移排列轨迹 ---
subplot(1, 2, 2);
hold on; grid on; box on;
title(['网格排列轨迹 (间距 ' num2str(gridSpacing) ' \mum)'], 'FontSize', 12, 'FontWeight', 'bold');
xlabel('X (\mum)', 'FontSize', 11);
ylabel('Y (\mum)', 'FontSize', 11);
axis equal;

% 4. 循环读取并绘制每个粒子的轨迹
for k = 1:numFiles
    filePath = fullfile(folderPath, csvFiles(k).name);
    
    % 读取 CSV 表格
    opts = detectImportOptions(filePath);
    opts.VariableNamingRule = 'preserve';
    tbl = readtable(filePath, opts);
    
    % 自动提取 X 和 Y 列
    colNames = tbl.Properties.VariableNames;
    xIdx = find(strcmpi(strtrim(colNames), 'X'), 1);
    yIdx = find(strcmpi(strtrim(colNames), 'Y'), 1);
    
    if isempty(xIdx) || isempty(yIdx)
        x_pix = tbl{:, 3};
        y_pix = tbl{:, 4};
    else
        x_pix = tbl{:, xIdx};
        y_pix = tbl{:, yIdx};
    end
    
    % 转换单位为微米 (\mum)
    x_um = x_pix * pix2um;
    y_um = y_pix * pix2um;
    
    % 计算瞬时速度 v (单位: \mum / s)
    dx = diff(x_um);
    dy = diff(y_um);
    dist = sqrt(dx.^2 + dy.^2); % 两帧间位移 (um)
    speed = dist / dt;          % 瞬时速度 (um / s)
    
    % 补齐数组长度与坐标点数一致
    speed_full = [speed; speed(end)];
    
    % 构造渐变线条画板数据 (插入 NaN 防止相邻轨迹首尾误连)
    x_patch_abs = [x_um; NaN];
    y_patch_abs = [y_um; NaN];
    v_patch     = [speed_full; NaN];
    
    % -------------------------------------------------------------
    % 绘制到子图 1：原始坐标（按速度着色）
    % -------------------------------------------------------------
    subplot(1, 2, 1);
    patch(x_patch_abs, y_patch_abs, v_patch, 'EdgeColor', 'interp', ...
        'FaceColor', 'none', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    % 起点 (黑色实心圆点) 与 终点 (白色填充方块)
    plot(x_um(1), y_um(1), 'o', 'MarkerSize', 3.5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot(x_um(end), y_um(end), 's', 'MarkerSize', 4.5, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
    
    % -------------------------------------------------------------
    % 绘制到子图 2：网格化平移（按速度着色）
    % -------------------------------------------------------------
    subplot(1, 2, 2);
    
    % 计算网格偏移量
    rowIdx = floor((k-1) / numCols);
    colIdx = mod(k-1, numCols);
    x_offset = colIdx * gridSpacing;
    y_offset = -rowIdx * gridSpacing;
    
    x_grid = (x_um - x_um(1)) + x_offset;
    y_grid = (y_um - y_um(1)) + y_offset;
    
    x_patch_grid = [x_grid; NaN];
    y_patch_grid = [y_grid; NaN];
    
    patch(x_patch_grid, y_patch_grid, v_patch, 'EdgeColor', 'interp', ...
        'FaceColor', 'none', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    plot(x_grid(1), y_grid(1), 'o', 'MarkerSize', 3.5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot(x_grid(end), y_grid(end), 's', 'MarkerSize', 4.5, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
    
    % 标注粒子编号
    if numFiles <= 60
        text(x_grid(1), y_grid(1) + gridSpacing*0.1, ['p=' num2str(k)], ...
            'Color', [0.2 0.2 0.2], 'FontSize', 8, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
end

% 5. 设置 Colorbar 刻度范围固定为 20 ~ 80 um/s
colormap(currentColormap);
c = colorbar;
c.Label.String = '瞬时速度 (\mum / s)';
c.Label.FontSize = 11;
c.Label.FontWeight = 'bold';

% 强制限制颜色区间为 [20, 80]
caxis([minSpeedLimit, maxSpeedLimit]); % 旧版本 MATLAB 适配
if exist('clim', 'file') || exist('clim', 'builtin')
    clim([minSpeedLimit, maxSpeedLimit]); % 新版本 MATLAB (R2022b+) 推荐语法
end

fprintf('绘制完成！Colorbar 色彩区间已固定为 %d ~ %d um/s。\n', minSpeedLimit, maxSpeedLimit);