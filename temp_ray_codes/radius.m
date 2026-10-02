% =========================================================================
% Ensemble Rotation Radius R(t) 0-4s (Y-Axis = 0-40 um)
% =========================================================================
clc; clear; close all;
warning('off', 'all');

%% 1. Parameters Setup
dataFolder = './';          % Folder path for CSV files
fps = 15;                   % Frame rate (frames/sec)
dt = 1 / fps;               % Time step (seconds)
scale = 0.48;               % Calibration ratio (um/px) - 1 px = 0.48 um
targetTime = 4.0;           % Target duration (seconds)
targetFrames = round(targetTime * fps) + 1; % 61 frames (0.0s to 4.0s)
nSelect = 20;               % Number of representative trajectories to select

fileList = dir(fullfile(dataFolder, '*.csv'));
numFiles = length(fileList);
if numFiles == 0
    error('No CSV trajectory files found! Please check the directory path.');
end
fprintf('Found %d trajectory files. Loading data...\n', numFiles);

%% 2. Data Extraction & Single-Particle Mean Radius Calculation
mat_R_all = NaN(targetFrames, numFiles);
traj_mean_R = NaN(numFiles, 1);

for k = 1:numFiles
    filePath = fullfile(fileList(k).folder, fileList(k).name);
    data = readtable(filePath);
    
    numFramesInFile = height(data);
    if numFramesInFile < 3
        continue;
    end
    
    currFrames = min(numFramesInFile, targetFrames);
    
    % Read coordinates and convert to microns
    if ismember('X', data.Properties.VariableNames)
        x = double(data.X(1:currFrames)) * scale;
        y = double(data.Y(1:currFrames)) * scale;
    else
        x = double(data{1:currFrames, 1}) * scale;
        y = double(data{1:currFrames, 2}) * scale;
    end
    
    % 3-point circumcircle fitting
    R = NaN(currFrames, 1);
    for i = 2:currFrames-1
        x1 = x(i-1); y1 = y(i-1);
        x2 = x(i);   y2 = y(i);
        x3 = x(i+1); y3 = y(i+1);
        
        a = hypot(x2 - x1, y2 - y1);
        b = hypot(x3 - x2, y3 - y2);
        c = hypot(x3 - x1, y3 - y1);
        area = 0.5 * abs(x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2));
        
        if area > 1e-5
            R(i) = (a * b * c) / (4 * area);
        end
    end
    
    % Edge padding and outlier filtering
    if currFrames >= 3
        R(1) = R(2);
        R(end) = R(end-1);
    end
    maxR_um = 150 * scale;
    R(R > maxR_um | R <= 0) = NaN;
    R = fillmissing(R, 'nearest');
    
    mat_R_all(1:currFrames, k) = R;
    traj_mean_R(k) = mean(R, 'omitnan');
end

%% 3. Select 20 Trajectories Closest to the Median
validIdx = find(~isnan(traj_mean_R));
validNum = length(validIdx);
if validNum == 0
    error('No valid trajectory data found.');
end

% Sort valid trajectories by mean radius
[~, sortOrder] = sort(traj_mean_R(validIdx));
sortedValidIdx = validIdx(sortOrder);

% Select nSelect (20) trajectories centered around the median
actualSelect = min(nSelect, validNum);
startRank = max(1, round(validNum/2 - actualSelect/2 + 1));
endRank = startRank + actualSelect - 1;
selectedIdx = sortedValidIdx(startRank:endRank);
mat_R = mat_R_all(:, selectedIdx);
fprintf('Selected %d representative trajectories around median index (out of %d total).\n', ...
    actualSelect, validNum);

%% 4. Calculate Ensemble Statistics on Selected Data
timeVec = (0:targetFrames-1)' * dt; % 0.0s to 4.0s
R_mean = mean(mat_R, 2, 'omitnan');
validCounts = sum(~isnan(mat_R), 2);
R_std  = std(mat_R, 0, 2, 'omitnan');
R_sem  = R_std ./ sqrt(validCounts);

%% 5. Plot Representative Ensemble Rotation Radius R(t)
figure('Color', 'w', 'Position', [100, 100, 900, 650]);
plotShadedError(timeVec, R_mean, R_sem, [0.12 0.47 0.71]);

% Title and Axis Labels with Extra Large Fonts
title('Ensemble Rotation Radius R(t) (Mean \pm SEM)', 'FontSize', 28, 'FontWeight', 'normal');
xlabel('Time t (s)', 'FontSize', 26, 'FontWeight', 'normal', 'Interpreter', 'tex');
ylabel('Radius R (\mum)', 'FontSize', 26, 'FontWeight', 'normal', 'Interpreter', 'tex');

% Fix Axis Limits: X = 0-4 s, Y = 0-40 um
xlim([0, 4]);
ylim([0, 40]);

% Axis Tick Font Size set to 24pt, LineWidth set to 1.0pt
set(gca, 'FontSize', 24, 'FontWeight', 'normal', 'LineWidth', 1.0);
grid off; box on;

%% Helper Function: Shaded Error Bar
function plotShadedError(x, meanVal, errVal, color)
    x = x(:)'; meanVal = meanVal(:)'; errVal = errVal(:)';
    upperBound = meanVal + errVal;
    lowerBound = meanVal - errVal;
    
    fill([x, fliplr(x)], [upperBound, fliplr(lowerBound)], color, ...
        'FaceAlpha', 0.25, 'EdgeColor', 'none', 'DisplayName', 'Mean \pm SEM');
    hold on;
    plot(x, meanVal, 'Color', color, 'LineWidth', 3.0, 'DisplayName', 'Mean R(t)');
    legend('Location', 'northeast', 'FontSize', 22);
end