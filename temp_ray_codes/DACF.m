% =========================================================================
% Directional Persistence & DACF Analysis Script (X = 0-4s, Y = -0.5-1)
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
dt = 1 / fps;               % Frame interval in seconds (1/15s)
scale = 0.48;               % Calibration ratio (um/px) - 1 pix = 0.48 um
targetTime = 4.0;           % Target max time delay (4 seconds)
maxTauLimit = round(targetTime * fps); % 60 frames (0 to 4 seconds)

fprintf('Processing %d CSV files for Directional Persistence...\n', numFiles);

%% 2. Calculate Directional Auto-Correlation Function (DACF)
mat_DACF = NaN(maxTauLimit, numFiles);
v0_array = zeros(numFiles, 1);

for k = 1:numFiles
    filePath = fullfile(fileList(k).folder, fileList(k).name);
    data = readtable(filePath);
    
    if ismember('X', data.Properties.VariableNames)
        x = data.X * scale; y = data.Y * scale;
    else
        x = data{:, 1} * scale; y = data{:, 2} * scale;
    end
    
    % Tangential displacements
    dx = diff(x); dy = diff(y);
    dr = sqrt(dx.^2 + dy.^2);
    
    v0_array(k) = mean(dr / dt, 'omitnan');
    
    % Unit direction vectors e(t) = [ex, ey]
    ex = dx ./ (dr + 1e-12);
    ey = dy ./ (dr + 1e-12);
    
    N = length(ex);
    tau_max_k = min(maxTauLimit, N - 1);
    
    for tau = 1:tau_max_k
        % Dot product: e(t) . e(t + tau)
        dot_prod = ex(1:end-tau) .* ex(1+tau:end) + ey(1:end-tau) .* ey(1+tau:end);
        mat_DACF(tau, k) = mean(dot_prod, 'omitnan');
    end
end

tauVec = (1:maxTauLimit)' * dt;
dacf_mean = mean(mat_DACF, 2, 'omitnan');
dacf_sem  = std(mat_DACF, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(mat_DACF), 2));

%% 3. Numerical Calculation for Persistence Metrics
% Keep metrics calculation for Subplot (b) without plotting fit curve
dacf_fit_fun = @(p, tau) cos(p(2) * tau) .* exp(-p(1) * tau);
p0 = [1.0, 5.0];

try
    p_fit = nlinfit(tauVec, dacf_mean, dacf_fit_fun, p0);
    Dr_fit = abs(p_fit(1));
    Omega0_fit = abs(p_fit(2));
catch
    p_exp = polyfit(tauVec(1:15), log(max(1e-3, dacf_mean(1:15))), 1);
    Dr_fit = -p_exp(1);
    Omega0_fit = 0;
end

tau_p = 1 / (Dr_fit + 1e-12);                 % Persistence time (s)
v0_mean = mean(v0_array);
L_p = v0_mean * tau_p;                        % Persistence length (um)

fprintf('\n===== DIRECTIONAL PERSISTENCE RESULTS =====\n');
fprintf('Rotational Diffusion Dr : %.4f s^-1\n', Dr_fit);
fprintf('Persistence Time tau_p  : %.4f s\n', tau_p);
fprintf('Persistence Length L_p  : %.4f um\n', L_p);

%% 4. Visualization
figure('Color', 'w', 'Position', [50, 80, 1250, 550]);

% --- Subplot (a): Ensemble DACF Decay Curve ---
subplot(1, 2, 1); hold on;
upperBand = dacf_mean + dacf_sem;
lowerBand = dacf_mean - dacf_sem;

% SEM Error Band
fill([tauVec; flipud(tauVec)], [upperBand; flipud(lowerBand)], ...
    [0.12 0.47 0.71], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'DisplayName', 'Mean \pm SEM');

% Ensemble Mean plotted as a smooth solid line
plot(tauVec, dacf_mean, 'Color', [0.12 0.47 0.71], 'LineWidth', 2.5, 'DisplayName', 'Ensemble DACF');

yline(0, 'k:', 'LineWidth', 1);
title('(a) Directional Auto-Correlation (DACF)', 'FontSize', 22, 'FontWeight', 'normal');
xlabel('Time Delay \tau (s)', 'FontSize', 20, 'FontWeight', 'normal', 'Interpreter', 'tex');
ylabel('Direction Correlation C(\tau)', 'FontSize', 20, 'FontWeight', 'normal', 'Interpreter', 'tex');

% --- Set Axis Limits: X = 0-4 s, Y = -0.5 to 1 ---
xlim([0, 4]);
ylim([-0.5, 1]);

set(gca, 'FontSize', 18, 'FontWeight', 'normal', 'LineWidth', 1.0);
grid off; box on;
legend('Location', 'northeast', 'FontSize', 14);

% --- Subplot (b): Persistence Length vs Speed ---
subplot(1, 2, 2); hold on;
scatter(v0_array, v0_array * tau_p, 60, [0.85 0.33 0.10], 'filled', 'MarkerAlpha', 0.7);
xline(v0_mean, 'r--', sprintf('Mean v_0 = %.2f', v0_mean), 'LineWidth', 1.5, 'FontSize', 14);

title('(b) Single Particle Persistence Length L_p', 'FontSize', 22, 'FontWeight', 'normal');
xlabel('Mean Speed v_0 (\mum / s)', 'FontSize', 20, 'FontWeight', 'normal', 'Interpreter', 'tex');
ylabel('Persistence Length L_p (\mum)', 'FontSize', 20, 'FontWeight', 'normal', 'Interpreter', 'tex');

set(gca, 'FontSize', 18, 'FontWeight', 'normal', 'LineWidth', 1.0);
grid off; box on;