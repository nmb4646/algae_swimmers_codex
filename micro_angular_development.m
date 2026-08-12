clear; clc; close all;

%%% SETTINGS
dir_w = './data/algae_motion/buffer_micro/';   % change if needed
dir_m = './data/algae_motion/mucus_micro/';    % change if needed

traj_ids = 1:20;

% File patterns
file_pattern_w = "Microalgae in water (%i).csv";   % change if needed
file_pattern_m = "Microalgae in mucus (%i).csv";    % change if needed

% Minimum allowed trajectory lengths (in frames)
min_len_w = 10;
min_len_m = 10;

% Number of normalized time points in [0,1]
nNorm = 200;

% Smoothing for theta if needed
use_smoothing = false;
sgolay_order = 3;
sgolay_window = 11;

% Plot style
col_w = [0.06 0.43 0.74];
col_m = [0.85 0.33 0.10];

%%% RUN BOTH CONDITIONS
[t_norm, theta_mean_w, theta_sem_w, n_used_w] = ...
    compute_ensemble_relative_theta(dir_w, file_pattern_w, traj_ids, min_len_w, nNorm, ...
    use_smoothing, sgolay_order, sgolay_window);

[~, theta_mean_m, theta_sem_m, n_used_m] = ...
    compute_ensemble_relative_theta(dir_m, file_pattern_m, traj_ids, min_len_m, nNorm, ...
    use_smoothing, sgolay_order, sgolay_window);

%%% PLOT
figure('Color','w'); hold on;

% Shaded SEM bands
fill([t_norm; flipud(t_norm)], ...
     [theta_mean_w-theta_sem_w; flipud(theta_mean_w+theta_sem_w)], ...
     col_w, 'FaceAlpha', 0.2, 'EdgeColor', 'none','HandleVisibility','off');

fill([t_norm; flipud(t_norm)], ...
     [theta_mean_m-theta_sem_m; flipud(theta_mean_m+theta_sem_m)], ...
     col_m, 'FaceAlpha', 0.2, 'EdgeColor', 'none','HandleVisibility','off');

% Means
plot(t_norm, theta_mean_w, '-', 'Color', col_w, 'LineWidth', 2.5, ...
    'DisplayName', sprintf('Water (n = %d)', n_used_w));
plot(t_norm, theta_mean_m, '-', 'Color', col_m, 'LineWidth', 2.5, ...
    'DisplayName', sprintf('Mucus (n = %d)', n_used_m));

grid on;
xlabel('Normalized time');
ylabel('\theta(t) - \theta(0)  [rad]');
title('Ensemble-averaged orientation, microalgae');
legend('Location','best');
set(gca,"FontSize",15)
set(gcf,"Position",[100,100,900,600])

%%% OPTIONAL: normalize each curve by its final value too
% Uncomment this block if you want orientation amplitude normalized as well
%{
[t_norm, theta_mean_w2, theta_sem_w2, n_used_w2] = ...
    compute_ensemble_relative_theta(dir_w, file_pattern_w, traj_ids, min_len_w, nNorm, ...
    use_smoothing, sgolay_order, sgolay_window, true);

[~, theta_mean_m2, theta_sem_m2, n_used_m2] = ...
    compute_ensemble_relative_theta(dir_m, file_pattern_m, traj_ids, min_len_m, nNorm, ...
    use_smoothing, sgolay_order, sgolay_window, true);

figure('Color','w'); hold on;
fill([t_norm; flipud(t_norm)], ...
     [theta_mean_w2-theta_sem_w2; flipud(theta_mean_w2+theta_sem_w2)], ...
     col_w, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([t_norm; flipud(t_norm)], ...
     [theta_mean_m2-theta_sem_m2; flipud(theta_mean_m2+theta_sem_m2)], ...
     col_m, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

plot(t_norm, theta_mean_w2, '-', 'Color', col_w, 'LineWidth', 2.5, ...
    'DisplayName', sprintf('Water (n = %d)', n_used_w2));
plot(t_norm, theta_mean_m2, '-', 'Color', col_m, 'LineWidth', 2.5, ...
    'DisplayName', sprintf('Mucus (n = %d)', n_used_m2));

grid on;
xlabel('Normalized time');
ylabel('Normalized relative orientation');
title('Ensemble-averaged normalized relative orientation');
legend('Location','best');
%}

%%% ---------- LOCAL FUNCTION ----------
function [t_norm, theta_mean, theta_sem, n_used] = ...
    compute_ensemble_relative_theta(dir_data, file_pattern, traj_ids, min_len, nNorm, ...
    use_smoothing, sgolay_order, sgolay_window, normalize_amplitude)

    if nargin < 9
        normalize_amplitude = false;
    end

    theta_mat = [];
    used_ids = [];

    warning('off','all');

    for n = traj_ids
        filename = dir_data + sprintf(file_pattern, n);

        if ~isfile(filename)
            fprintf('Skipping missing file: %s\n', filename);
            continue;
        end

        T = readtable(filename);

        if ~all(ismember({'X','Y'}, T.Properties.VariableNames))
            fprintf('Skipping file without X/Y columns: %s\n', filename);
            continue;
        end

        x = T.X(:);
        y = T.Y(:);

        if numel(x) < min_len
            fprintf('Skipping traj %d: length %d < min_len %d\n', n, numel(x), min_len);
            continue;
        end

        % Orientation from trajectory tangent
        dx = gradient(x);
        dy = gradient(y);
        theta = atan2(dy, dx);
        theta_unwrapped = unwrap(theta);

        if use_smoothing && numel(theta_unwrapped) >= sgolay_window
            theta_unwrapped = sgolayfilt(theta_unwrapped, sgolay_order, sgolay_window);
        end

        % Relative orientation
        theta_rel = theta_unwrapped - theta_unwrapped(1);

        % Normalize time to [0,1]
        t_this = linspace(0, 1, numel(theta_rel))';
        t_norm = linspace(0, 1, nNorm)';

        % Interpolate onto common normalized time grid
        theta_interp = interp1(t_this, theta_rel, t_norm, 'linear');

        % Optional amplitude normalization
        if normalize_amplitude
            denom = theta_interp(end);
            if abs(denom) > eps
                theta_interp = theta_interp / denom;
            else
                continue;
            end
        end

        theta_mat = [theta_mat, theta_interp];
        used_ids(end+1) = n; %#ok<AGROW>
    end

    warning('on','all');

    if isempty(theta_mat)
        error('No valid trajectories found in %s', dir_data);
    end

    theta_mean = mean(theta_mat, 2, 'omitnan');
    theta_sem  = std(theta_mat, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(theta_mat), 2));
    n_used = size(theta_mat, 2);

    fprintf('\nDataset: %s\n', dir_data);
    fprintf('Used trajectories: %s\n', mat2str(used_ids));
    fprintf('n used: %d\n', n_used);
end