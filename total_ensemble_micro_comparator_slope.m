%clear; clc; close all;

%%%% SETTINGS
dir_w = './data_standardized/micro_water/';
dir_m = './data_standardized/micro_mucus/';

use_mucus = true;                  % true = mucus, false = water
if use_mucus
    traj_ids = 1:51;
else
    traj_ids = 1:50;
end

use_smoothing = true;
go = 3;
gw = 11;
common_dt = [];                    % [] = use coarsest native dt among loaded trajectories
resample_method = 'pchip';

alpha_win = 5;                     % sliding window size on log-log axes
min_pairs = 3;                     % minimum pair count required for slope estimate
frac_latefit = [0.50 0.70];        % late window for a single summary alpha
min_fit_pts = 3;

%%%% CHOOSE DATASET
if use_mucus
    dir_data = dir_m;
    file_pattern = "Microalgae in mucus (%i).csv";
    cond_name = "mucus";
else
    dir_data = dir_w;
    file_pattern = "Microalgae in water (%i).csv";
    cond_name = "water";
end

[traj_data, track_lengths, used_ids, native_dts, common_dt] = ...
    load_standardized_resampled_trajectories(dir_data, file_pattern, traj_ids, ...
    use_smoothing, go, gw, common_dt, resample_method);

nTraj = numel(traj_data);

fprintf('Condition: %s\n', cond_name);
fprintf('Loaded %d trajectories\n', nTraj);
fprintf('Trajectory IDs: %s\n', mat2str(used_ids'));
fprintf('Common resampled dt: %.6g s\n', common_dt);
fprintf('Native dt range: [%.6g, %.6g] s\n', min(native_dts), max(native_dts));

%%%% ENSEMBLE MSD UP TO LONGEST AVAILABLE LAG
maxLag = max(track_lengths) - 1;

tau = (1:maxLag)' * common_dt;
msd_sum = zeros(maxLag,1);
pair_count = zeros(maxLag,1);

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
    end
end

ensemble_msd = msd_sum ./ pair_count;
ensemble_msd(pair_count == 0) = NaN;

%%%% LOCAL BEST-FIT SLOPE ON LOG-LOG AXES
valid = isfinite(ensemble_msd) & ensemble_msd > 0 & pair_count >= min_pairs;
tau_valid = tau(valid);
msd_valid = ensemble_msd(valid);
nValid = numel(tau_valid);

if nValid < min_fit_pts
    error('Not enough valid MSD points to estimate local slopes.');
end

logtau = log10(tau_valid);
logmsd = log10(msd_valid);

if mod(alpha_win, 2) == 0
    error('alpha_win must be odd.');
end

halfwin = floor(alpha_win / 2);
alpha_tau = nan(size(tau_valid));

for i = 1:nValid
    i1 = max(1, i - halfwin);
    i2 = min(nValid, i + halfwin);

    if (i2 - i1 + 1) >= min_fit_pts
        p = polyfit(logtau(i1:i2), logmsd(i1:i2), 1);
        alpha_tau(i) = p(1);
    end
end

%%%% LATE-TIME SUMMARY SLOPE
to_idx = @(fr) max(1, min(nValid, round(fr * nValid)));
idx_late = to_idx(frac_latefit(1)) + 1 : to_idx(frac_latefit(2));

if numel(idx_late) < min_fit_pts
    error('Late-fit window contains too few points. Adjust frac_latefit.');
end

p_late = polyfit(logtau(idx_late), logmsd(idx_late), 1);
alpha_late = p_late(1);

fprintf('Late-window fitted alpha = %.3f\n', alpha_late);

%%%% ERROR ESTIMATE
rel_err = nan(maxLag,1);
Meff = nan(maxLag,1);

for lag = 1:maxLag
    idx = track_lengths > lag;
    Meff(lag) = sum(1.5 * track_lengths(idx) / lag);
    rel_err(lag) = 1 / sqrt(Meff(lag));   % for 2D diffusion
end

percent_err = 100 * rel_err;

%%%% PLOT
ax1 = subplot(2,1,1);
hold on;
grid on;
semilogx(tau_valid, alpha_tau, 'o-', 'LineWidth', 2, 'DisplayName', cond_name + ", micro");
yline(2, '--', 'Ballistic', 'HandleVisibility', 'off');
yline(1, '--', 'Diffusive', 'HandleVisibility', 'off');
yline(0, '--', 'Confined', 'HandleVisibility', 'off');
yline(alpha_late, ':', sprintf('Late fit = %.2f', alpha_late), ...
    'LineWidth', 1.5, 'DisplayName', sprintf('Late fit = %.2f', alpha_late));
xlabel('\tau (s)');
ylabel('Best-fit slope \alpha(\tau)');
title('Ensemble MSD Slope Comparison, microalgae');
legend('Location', 'best');

ax2 = subplot(2,1,2);
hold on;
grid on;
plot(tau, percent_err, '-', 'LineWidth', 2.5);
ylabel('Error %');
xlabel('\tau (s)');
title('Error quantification');
set(gca, 'XScale', 'log');

linkaxes([ax1, ax2], 'x');
xlim(ax1, [min(tau_valid) max(tau_valid)]);

set(gcf, 'Position', [100,100,1400,900]);
subplot(2,1,1); set(gca, 'FontSize', 18);
subplot(2,1,2); set(gca, 'FontSize', 18);

%%%% SUMMARY
fprintf('Shortest trajectory: %d frames\n', min(track_lengths));
fprintf('Longest trajectory: %d frames\n', max(track_lengths));
fprintf('Max lag analyzed: %d resampled steps = %.3f s\n', maxLag, maxLag * common_dt);
fprintf('Valid slope points: %d / %d\n', nValid, numel(tau));
fprintf('Sliding alpha window: %d points\n', alpha_win);
fprintf('Late fit window: %d points\n', numel(idx_late));
