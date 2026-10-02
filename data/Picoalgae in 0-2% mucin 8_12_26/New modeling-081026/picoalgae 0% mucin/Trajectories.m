% =========================================================
% MATLAB 轨迹绘制 (防错修复版)
% 参数: Pix to Micron = 0.48, FPS = 15, Speed Colorbar = 20-80 um/s
% =========================================================

clear; clc; close all;

% 1. 参数定义
pix2um = 0.48;       % 像素转微米系数 (um / pixel)
fps = 15;            % 帧率 (fps)
dt = 1 / fps;        % 帧时间间隔 (s)
folderPath = './';   % CSV 文件所在文件夹
gridSpacing = 150;   % 网格排列间距 (um)

minSpeedLimit = 20;  % Colorbar 速度下限 (um/s)
maxSpeedLimit = 80;  % Colorbar 速度上限 (um/s)

% 2. 检索并严格筛选纯数字命名的 CSV 文件 (避免误读 summary.csv 等文件)
filePattern = fullfile(folderPath, '*.csv');
allFiles = dir(filePattern);

if isempty(allFiles)
    error('未在指定路径下找到任何 CSV 文件，请检查文件夹路径！');
end

% 过滤：只保留文件名中包含有效阿拉伯数字的文件
validFiles = {};
fileNums = [];

for i = 1:length(allFiles)
    % 提取文件名中的数字
    numMatch = regexp(allFiles(i).name, '^\d+(?=\.csv$)', 'match'); 
    if ~isempty(numMatch)
        validFiles{end+1} = allFiles(i).name; %#ok<SAGROW>
        fileNums(end+1) = str2double(numMatch{1}); %#ok<SAGROW>
    end
end

if isempty(validFiles)
    error('未找到纯数字命名的 CSV 文件（如 1.csv, 2.csv），请检查文件名！');
end

% 按数值升序排列
[~, sortIdx] = sort(fileNums);
csvFiles = validFiles(sortIdx);
numFiles = length(csvFiles);

fprintf('成功找到 %d 个轨迹文件，开始处理...\n', numFiles);

% 3. 计算网格行列数
numCols = ceil(sqrt(numFiles)); 
numRows = ceil(numFiles / numCols); 

% 4. 初始化窗口
fig = figure('Name', '轨迹绘制 (防错修复版)', 'Color', 'w', 'Position', [100, 100, 1250, 550]);

ax1 = subplot(1, 2, 1);
hold(ax1, 'on'); grid(ax1, 'on'); box(ax1, 'on');
title(ax1, '真实物理位置轨迹 (颜色表示瞬时速度)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel(ax1, 'X (\mum)', 'FontSize', 11);
ylabel(ax1, 'Y (\mum)', 'FontSize', 11);
axis(ax1, 'equal');

ax2 = subplot(1, 2, 2);
hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');
title(ax2, ['网格排列轨迹 (间距 ' num2str(gridSpacing) ' \mum)'], 'FontSize', 12, 'FontWeight', 'bold');
xlabel(ax2, 'X (\mum)', 'FontSize', 11);
ylabel(ax2, 'Y (\mum)', 'FontSize', 11);
axis(ax2, 'equal');

% 5. 循环读取并绘制
for k = 1:numFiles
    fileName = csvFiles{k};
    filePath = fullfile(folderPath, fileName);
    
    try
        % 读取 CSV 文件
        opts = detectImportOptions(filePath);
        opts.VariableNamingRule = 'preserve';
        tbl = readtable(filePath, opts);
        
        % 获取表格列名
        colNames = tbl.Properties.VariableNames;
        
        % 智能匹配 X 和 Y 列名
        xIdx = find(~cellfun(@isempty, regexpi(strtrim(colNames), '^x$|^x\b')), 1);
        yIdx = find(~cellfun(@isempty, regexpi(strtrim(colNames), '^y$|^y\b')), 1);
        
        % 如果智能匹配失败，回退至第 3 列 (X) 和 第 4 列 (Y)
        if isempty(xIdx) || isempty(yIdx)
            if width(tbl) >= 4
                xIdx = 3;
                yIdx = 4;
            else
                warning('文件 %s 列数不满足要求(当前仅 %d 列)，已跳过该文件。', fileName, width(tbl));
                continue; % 安全跳过该文件
            end
        end
        
        % 提取坐标数据
        x_pix = tbl{:, xIdx};
        y_pix = tbl{:, yIdx};
        
        % 转换为微米 (\mum)
        x_um = x_pix * pix2um;
        y_um = y_pix * pix2um;
        
        % 计算瞬时速度 v (um / s)
        dx = diff(x_um);
        dy = diff(y_um);
        dist = sqrt(dx.^2 + dy.^2); 
        speed = dist / dt;          
        speed_full = [speed; speed(end)]; % 补齐维度
        
        % 构造 Patch 数据
        x_patch_abs = [x_um; NaN];
        y_patch_abs = [y_um; NaN];
        v_patch     = [speed_full; NaN];
        
        % --- 绘制子图 1 (真实空间) ---
        patch('Parent', ax1, 'XData', x_patch_abs, 'YData', y_patch_abs, 'CData', v_patch, ...
            'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 1.2);
        plot(ax1, x_um(1), y_um(1), 'o', 'MarkerSize', 3.5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
        plot(ax1, x_um(end), y_um(end), 's', 'MarkerSize', 4.5, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
        
        % --- 绘制子图 2 (网格排列) ---
        rowIdx = floor((k-1) / numCols);
        colIdx = mod(k-1, numCols);
        x_offset = colIdx * gridSpacing;
        y_offset = -rowIdx * gridSpacing;
        
        x_grid = (x_um - x_um(1)) + x_offset;
        y_grid = (y_um - y_um(1)) + y_offset;
        
        x_patch_grid = [x_grid; NaN];
        y_patch_grid = [y_grid; NaN];
        
        patch('Parent', ax2, 'XData', x_patch_grid, 'YData', y_patch_grid, 'CData', v_patch, ...
            'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 1.2);
        plot(ax2, x_grid(1), y_grid(1), 'o', 'MarkerSize', 3.5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
        plot(ax2, x_grid(end), y_grid(end), 's', 'MarkerSize', 4.5, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
        
        % 标注文件编号
        if numFiles <= 60
            text(ax2, x_grid(1), y_grid(1) + gridSpacing*0.1, ['p=' num2str(fileNums(k))], ...
                'Color', [0.2 0.2 0.2], 'FontSize', 8, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        end
        
    catch ME
        warning('读取文件 %s 时发生错误: %s，已自动跳过。', fileName, ME.message);
    end
end

% 6. 设置色彩与 Colorbar
colormap(fig, turbo);

caxis(ax1, [minSpeedLimit, maxSpeedLimit]);
caxis(ax2, [minSpeedLimit, maxSpeedLimit]);

c = colorbar(ax2, 'eastoutside');
c.Label.String = '瞬时速度 (\mum / s)';
c.Label.FontSize = 11;
c.Label.FontWeight = 'bold';

fprintf('完成！已顺利绘制全部有效轨迹。\n');