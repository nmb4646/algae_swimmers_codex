clear; clc; close all;

%%%% SETTINGS
base_dir = './data/Picoalgae in 0-2% mucin/';

condition_dirs = [ ...
    "picoalgae 0% mucin"; ...
    "picoalgae 1% mucin"; ...
    "picoalgae 2% mucin"];

condition_names = ["0% mucin", "1% mucin", "2% mucin"];

fps = 15;
native_dt = 1 / fps;
pix_m_per_px = 0.48e-6;

use_smoothing = true;
go = 3;
gw = 11;
msd_mean_type = 'geometric';      % 'arithmetic' or 'geometric'
velo_mean_type = 'geometric';     % 'arithmetic' or 'geometric'
use_resampling = false;            % true = pchip-resample FrameNumber gaps to uniform native_dt
show_late_fit = false;             % dashed overlay for free-alpha late fit
show_brownian_line = false;        % plot theoretical MSD = D0*tau line

% -------- tunable fit windows (fractions of valid MSD points) --------
frac_ballistic = [0.00 0.03];      % alpha = 2 fit window
frac_diffusive = [0.03 0.15];      % alpha = 1 fit window
frac_latefit   = [0.40 1];      % free alpha fit window
min_fit_pts    = 3;                % minimum points required in a fit window
% --------------------------------------------------------------------

colors = lines(numel(condition_names));

figure('Color','w');
ax1 = subplot(2,1,1); hold on; grid on;
ax2 = subplot(2,1,2); hold on; grid on;

all_tau_valid = [];

for c = 1:numel(condition_names)
    dir_data = fullfile(base_dir, condition_dirs(c));
    cond_name = condition_names(c);

    [traj_data, track_lengths, used_files] = load_raw_pico_mucin_trajectories( ...
        dir_data, pix_m_per_px, native_dt, use_smoothing, go, gw, use_resampling);

    nTraj = numel(traj_data);

    fprintf('\nCondition: %s\n', cond_name);
    fprintf('MSD mean type: %s\n', msd_mean_type);
    fprintf('Velocity mean type: %s\n', velo_mean_type);
    fprintf('Resampling enabled: %d\n', use_resampling);
    fprintf('Loaded %d trajectories\n', nTraj);
    fprintf('Files: %s\n', strjoin(used_files, ", "));
    fprintf('Native dt: %.6g s\n', native_dt);

    [mean_speeds, angular_velocities] = compute_pico_transport_metrics(traj_data, native_dt);
    mean_speed = aggregate_positive(mean_speeds, velo_mean_type);
    mean_angular_velocity = aggregate_positive(abs(angular_velocities), velo_mean_type);
    d_eff = D2(1e-6, mean_speed, mean_angular_velocity);

    fprintf('%s mean speed = %.3e m/s\n', titlecase_mean_type(velo_mean_type), mean_speed);
    fprintf('%s mean |angular velocity| = %.3e rad/s\n', titlecase_mean_type(velo_mean_type), mean_angular_velocity);
    fprintf('D2 effective diffusion from %s mean speed/omega = %.3e m^2/s\n', velo_mean_type, d_eff);

    %%%% ENSEMBLE MSD UP TO LONGEST AVAILABLE LAG
    maxLag = max(track_lengths) - 1;

    tau = (1:maxLag)' * native_dt;
    msd_sum = zeros(maxLag,1);
    pair_count = zeros(maxLag,1);
    log_msd_sum = zeros(maxLag,1);
    log_pair_count = zeros(maxLag,1);

    for k = 1:nTraj
        xy = traj_data{k};
        x = xy(:,1);
        y = xy(:,2);
        N = numel(x);

        for lag = 1:(N-1)
            dx = x(1+lag:N) - x(1:N-lag);
            dy = y(1+lag:N) - y(1:N-lag);
            sqdisp = dx.^2 + dy.^2;

            msd_sum(lag) = msd_sum(lag) + sum(sqdisp);
            pair_count(lag) = pair_count(lag) + numel(sqdisp);

            positive_sqdisp = sqdisp(sqdisp > 0);
            log_msd_sum(lag) = log_msd_sum(lag) + sum(log(positive_sqdisp));
            log_pair_count(lag) = log_pair_count(lag) + numel(positive_sqdisp);
        end
    end

    switch lower(msd_mean_type)
        case 'arithmetic'
            ensemble_msd = msd_sum ./ pair_count;
            ensemble_msd(pair_count == 0) = NaN;
        case 'geometric'
            ensemble_msd = exp(log_msd_sum ./ log_pair_count);
            ensemble_msd(log_pair_count == 0) = NaN;
        otherwise
            error('Unsupported msd_mean_type: %s', msd_mean_type);
    end

    valid = isfinite(ensemble_msd) & ensemble_msd > 0 & pair_count > 0;
    tau_valid = tau(valid);
    msd_valid = ensemble_msd(valid);
    nValid = numel(tau_valid);

    if nValid < min_fit_pts
        error('Not enough valid MSD points to fit for %s.', cond_name);
    end

    to_idx = @(fr) max(1, min(nValid, round(fr * nValid)));

    idx_ball = to_idx(frac_ballistic(1)) + 1 : to_idx(frac_ballistic(2));
    idx_diff = to_idx(frac_diffusive(1)) + 1 : to_idx(frac_diffusive(2));
    idx_late = to_idx(frac_latefit(1)) + 1 : to_idx(frac_latefit(2));

    if numel(idx_ball) < min_fit_pts || numel(idx_diff) < min_fit_pts || numel(idx_late) < min_fit_pts
        error('One or more fit windows contain too few points for %s. Adjust the fraction settings.', cond_name);
    end

    logtau = log10(tau_valid);
    logmsd = log10(msd_valid);

    b_ball = mean(logmsd(idx_ball) - 2*logtau(idx_ball));
    fit_ball = 10.^(2*logtau(idx_ball) + b_ball);

    b_diff = mean(logmsd(idx_diff) - logtau(idx_diff));
    fit_diff = 10.^(logtau(idx_diff) + b_diff);

    p_late = polyfit(logtau(idx_late), logmsd(idx_late), 1);
    alpha_late = p_late(1);
    K_late = 10^(p_late(2));
    fit_late = 10.^(polyval(p_late, logtau(idx_late)));

    fprintf('Late-window fitted alpha = %.3f\n', alpha_late);
    fprintf('Late-window fitted K = %.3e m^2 / s^alpha\n', K_late);

    loglog(ax1, tau, ensemble_msd, 'o-', ...
        'LineWidth', 2, ...
        'Color', colors(c,:), ...
        'DisplayName', cond_name);

    if show_late_fit
        loglog(ax1, tau_valid(idx_late), fit_late, '--', ...
            'LineWidth', 2.5, ...
            'Color', colors(c,:), ...
            'DisplayName', sprintf('%s late fit, \\alpha=%.2f', cond_name, alpha_late));
    end

    rel_err = nan(maxLag,1);
    Meff = nan(maxLag,1);

    for lag = 1:maxLag
        idx = track_lengths > lag;
        Meff(lag) = sum(1.5 * track_lengths(idx) / lag);
        rel_err(lag) = 1 / sqrt(Meff(lag));   % for 2D diffusion
    end

    percent_err = 100 * rel_err;

    plot(ax2, tau, percent_err, '-', ...
        'LineWidth', 2.5, ...
        'Color', colors(c,:), ...
        'DisplayName', cond_name);

    all_tau_valid = [all_tau_valid; tau_valid(:)]; %#ok<AGROW>

    fprintf('Shortest trajectory: %d frames\n', min(track_lengths));
    fprintf('Longest trajectory: %d frames\n', max(track_lengths));
    fprintf('Max lag analyzed: %d steps = %.3f s\n', maxLag, maxLag*native_dt);
    fprintf('Ballistic fit window: %d points\n', numel(idx_ball));
    fprintf('Diffusive fit window: %d points\n', numel(idx_diff));
    fprintf('Late fit window: %d points\n', numel(idx_late));
end

xlabel(ax1, '\tau (s)');
ylabel(ax1, 'Ensemble MSD(\tau) (m^2)');
title(ax1, 'Ensemble MSD Comparison, picoalgae in 0-2% mucin');
legend(ax1, 'Location', 'best');
set(ax1, 'XScale', 'log');
set(ax1, 'YScale', 'log');

xlabel(ax2, '\tau (s)');
ylabel(ax2, 'Error %');
title(ax2, 'Error quantification');
set(ax2, 'XScale', 'log');
legend(ax2, 'Location', 'best');

linkaxes([ax1, ax2], 'x');
if ~isempty(all_tau_valid)
    xlim(ax1, [min(all_tau_valid), max(all_tau_valid)]);

    if show_brownian_line
        D0 = brownian_diffusion_coefficient(1e-6);
        tau_ref = [min(all_tau_valid), max(all_tau_valid)];
        loglog(ax1, tau_ref, D0 * tau_ref, 'k:', ...
            'LineWidth', 2, ...
            'DisplayName', sprintf('MSD = D_0\\tau, D_0=%.2e m^2/s', D0));
        legend(ax1, 'Location', 'best');
    end
end
ylim(ax2, [1, 50]);

set(gcf, 'Position', [100, 100, 1400, 900]);
set(ax1, 'FontSize', 15);
set(ax2, 'FontSize', 15);

function [traj_data, track_lengths, used_files] = load_raw_pico_mucin_trajectories( ...
    dir_data, pix_m_per_px, native_dt, use_smoothing, go, gw, use_resampling)

    if ~isfolder(dir_data)
        error('Input folder not found: %s', dir_data);
    end

    files = dir(fullfile(dir_data, '*.csv'));
    files = files(~[files.isdir]);
    if isempty(files)
        error('No CSV files found in %s.', dir_data);
    end

    file_nums = nan(numel(files), 1);
    for k = 1:numel(files)
        [~, name, ~] = fileparts(files(k).name);
        file_nums(k) = str2double(name);
    end
    [~, order] = sort(file_nums);
    files = files(order);

    traj_data = {};
    track_lengths = [];
    used_files = strings(0, 1);

    for k = 1:numel(files)
        filename = fullfile(files(k).folder, files(k).name);

        warning('off', 'all');
        T = readtable(filename);
        warning('on', 'all');

        vars = string(T.Properties.VariableNames);
        if ~all(ismember(["X", "Y"], vars))
            fprintf('Missing X/Y columns, skipping: %s\n', filename);
            continue;
        end

        x_px = T.X(:);
        y_px = T.Y(:);

        if ismember("FrameNumber", vars)
            frame_number = T.FrameNumber(:);
        else
            frame_number = (0:numel(x_px)-1)';
        end

        valid_rows = isfinite(frame_number) & isfinite(x_px) & isfinite(y_px);
        x_px = x_px(valid_rows);
        y_px = y_px(valid_rows);
        frame_number = frame_number(valid_rows);

        if numel(x_px) < 2
            fprintf('Trajectory too short, skipping: %s\n', filename);
            continue;
        end

        [frame_number, unique_idx] = unique(frame_number, 'stable');
        x_px = x_px(unique_idx);
        y_px = y_px(unique_idx);

        x = (x_px - x_px(1)) * pix_m_per_px;
        y = (y_px - y_px(1)) * pix_m_per_px;
        t_native = (frame_number - frame_number(1)) * native_dt;

        if use_smoothing && numel(x) >= gw
            x = sgolayfilt(x, go, gw);
            y = sgolayfilt(y, go, gw);
        end

        if any(~isfinite(t_native)) || any(~isfinite(x)) || any(~isfinite(y))
            fprintf('Non-finite values after preprocessing, skipping: %s\n', filename);
            continue;
        end

        if use_resampling
            duration = t_native(end) - t_native(1);
            n_uniform = floor(duration / native_dt) + 1;
            if n_uniform < 2
                fprintf('Trajectory too short after resampling, skipping: %s\n', filename);
                continue;
            end

            t_uniform = t_native(1) + (0:n_uniform-1)' * native_dt;
            x = interp1(t_native, x, t_uniform, 'pchip');
            y = interp1(t_native, y, t_uniform, 'pchip');

            if any(~isfinite(x)) || any(~isfinite(y))
                fprintf('Interpolation produced non-finite values, skipping: %s\n', filename);
                continue;
            end
        end

        traj_data{end+1,1} = [x(:), y(:)]; %#ok<AGROW>
        track_lengths(end+1,1) = numel(x); %#ok<AGROW>
        used_files(end+1,1) = string(files(k).name); %#ok<AGROW>
    end

    if isempty(traj_data)
        error('No trajectories were loaded from %s.', dir_data);
    end
end

function [mean_speeds, angular_velocities] = compute_pico_transport_metrics(traj_data, dt)
    nTraj = numel(traj_data);
    mean_speeds = nan(nTraj, 1);
    angular_velocities = nan(nTraj, 1);

    for k = 1:nTraj
        xy = traj_data{k};
        x = xy(:,1);
        y = xy(:,2);

        if numel(x) < 2
            continue;
        end

        dx = gradient(x) / dt;
        dy = gradient(y) / dt;
        speed = hypot(dx, dy);

        theta = atan2(dy, dx);
        theta_unwrapped = unwrap(theta);
        total_turn = theta_unwrapped(end) - theta_unwrapped(1);

        mean_speeds(k) = mean(speed, 'omitnan');
        angular_velocities(k) = total_turn / (numel(x) * dt);
    end
end

function value = aggregate_positive(values, mean_type)
    values = values(:);
    values = values(isfinite(values) & values > 0);

    if isempty(values)
        value = NaN;
        return;
    end

    switch lower(string(mean_type))
        case "arithmetic"
            value = mean(values);
        case "geometric"
            value = exp(mean(log(values)));
        otherwise
            error('Unsupported velo_mean_type: %s', mean_type);
    end
end

function label = titlecase_mean_type(mean_type)
    switch lower(string(mean_type))
        case "arithmetic"
            label = "Arithmetic";
        case "geometric"
            label = "Geometric";
        otherwise
            error('Unsupported velo_mean_type: %s', mean_type);
    end
end

function D0 = brownian_diffusion_coefficient(a)
    kb = 1.38e-23;
    T = 293;
    mu = 0.001;

    D0 = kb * T / (6 * pi * mu * a);
end
