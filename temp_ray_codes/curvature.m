% =========================================================================
% Ensemble Trajectory Curvature Statistical Analysis (Y-Axis = 0-12)
% =========================================================================
clc; clear; close all;
warning('off', 'MATLAB:table:ModifiedAndSavedVariableNames');

%% 1. Parameters Setup
dataFolder = './';          % Path to CSV files
fileList = dir(fullfile(dataFolder, '*.csv'));
numFiles = length(fileList);
if numFiles == 0
    error('No CSV files found in directory!');
end

% --- Calibration Parameters ---
fps = 15;                   % Frame rate (frames/sec)
dt = 1 / fps;               % Frame interval (s)
scale = 0.48;               % Calibration ratio (um/px) - 1 px = 0.48 um

% Grid for curvature density alignment (in um^-1)
max_k_grid = 0.5;           % Max curvature limit (1/um)
grid_k = linspace(0, max_k_grid, 200); 
pdf_matrix = zeros(numFiles, length(grid_k)); 
particle_mean_k   = zeros(numFiles, 1); 
particle_median_k = zeros(numFiles, 1); 

fprintf('Processing %d CSV trajectory files...\n', numFiles);

%% 2. Batch Processing & Density Estimation
for k = 1:numFiles
    filePath = fullfile(fileList(k).folder, fileList(k).name);
    data = readtable(filePath);
    
    if ismember('X', data.Properties.VariableNames)
        x = double(data.X) * scale; 
        y = double(data.Y) * scale;
    else
        x = double(data{:, 1}) * scale; 
        y = double(data{:, 2}) * scale;
    end
    
    N = length(x);
    
    % Three-point Menger curvature calculation (in um^-1)
    k_val = zeros(N, 1);
    for i = 2:N-1
        a = hypot(x(i)-x(i-1), y(i)-y(i-1));
        b = hypot(x(i+1)-x(i), y(i+1)-y(i));
        c = hypot(x(i+1)-x(i-1), y(i+1)-y(i-1));
        area = 0.5 * abs(x(i-1)*(y(i)-y(i+1)) + x(i)*(y(i+1)-y(i-1)) + x(i+1)*(y(i-1)-y(i)));
        if (a * b * c) > 1e-6
            k_val(i) = (4 * area) / (a * b * c);
        end
    end
    
    % Filter invalid values and outliers
    valid_k = k_val(k_val > 0 & k_val < max_k_grid);
    
    if ~isempty(valid_k)
        pdf_matrix(k, :) = ksdensity(valid_k, grid_k);
        particle_mean_k(k)   = mean(valid_k);
        particle_median_k(k) = median(valid_k);
    end
end

%% 3. Compute Ensemble Statistics
pdf_mean = mean(pdf_matrix, 1, 'omitnan');
pdf_std  = std(pdf_matrix, 0, 1, 'omitnan');
pdf_sem  = pdf_std / sqrt(numFiles);

%% 4. Visualization
figure('Color', 'w', 'Position', [50, 80, 1150, 480]);

% --- Subplot (a): Ensemble Curvature PDF ---
subplot(1, 2, 1); hold on;
upper_sem = pdf_mean + pdf_sem;
lower_sem = max(0, pdf_mean - pdf_sem);

% SEM Error Band
fill([grid_k, fliplr(grid_k)], [upper_sem, fliplr(lower_sem)], ...
    [0.12 0.47 0.71], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Mean \pm SEM');

% Ensemble Mean PDF Line
plot(grid_k, pdf_mean, 'Color', [0.12 0.47 0.71], 'LineWidth', 2.5, 'DisplayName', 'Ensemble Mean PDF');

% Dominant curvature vertical line (Clean without text)
[max_y, max_idx] = max(pdf_mean);
peak_k = grid_k(max_idx);
xline(peak_k, 'r--', 'LineWidth', 1.2);

title(sprintf('(a) Curvature Density PDF (N = %d)', numFiles), 'FontSize', 20, 'FontWeight', 'normal');
xlabel('Curvature \kappa (\mum^{-1})', 'FontSize', 18, 'FontWeight', 'normal', 'Interpreter', 'tex');
ylabel('Probability Density', 'FontSize', 18, 'FontWeight', 'normal', 'Interpreter', 'tex');

% --- Axis Limits: X = 0-0.3, Y = 0-12 ---
xlim([0, 0.3]);
ylim([0, 12]);

set(gca, 'FontSize', 16, 'LineWidth', 1.0);
legend('Location', 'northeast', 'FontSize', 12);
grid off; box on;

% --- Subplot (b): Individual Variance (Scatter + Boxplot + Error Bar) ---
subplot(1, 2, 2); hold on;

% 1. Scatter with horizontal jitter
jitter = (rand(numFiles, 1) - 0.5) * 0.15;
scatter(1 + jitter, particle_mean_k, 30, [0.6 0.6 0.6], 'filled', 'MarkerAlpha', 0.6);

% 2. Boxplot
hBox = boxchart(ones(numFiles, 1), particle_mean_k, 'BoxWidth', 0.3, 'CapStyle', 'none');
hBox.BoxFaceColor = [0.85 0.33 0.10];
hBox.BoxFaceAlpha = 0.4;
hBox.MarkerStyle = 'none';

% 3. Mean +- STD Error Bar
pop_mean = mean(particle_mean_k);
pop_std  = std(particle_mean_k);
errorbar(1.3, pop_mean, pop_std, 'k', 'LineWidth', 1.5, 'CapSize', 8, 'Marker', 'o', ...
    'MarkerFaceColor', 'r', 'MarkerSize', 6);

set(gca, 'XTick', 1, 'XTickLabel', {'Particle Mean \kappa'}, 'FontSize', 16, 'LineWidth', 1.0);
title('(b) Individual Curvature Variance', 'FontSize', 20, 'FontWeight', 'normal');
ylabel('Mean Curvature \kappa (\mum^{-1})', 'FontSize', 18, 'FontWeight', 'normal', 'Interpreter', 'tex');
xlim([0.6, 1.7]);
grid off; box on;