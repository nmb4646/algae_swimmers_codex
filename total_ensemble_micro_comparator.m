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
msd_mean_type = 'geometric';      % 'arithmetic' or 'geometric'

% -------- tunable fit windows (fractions of valid MSD points) --------
frac_ballistic = [0 0.10];      % alpha = 2 fit window
frac_diffusive = [0.10 0.50];      % alpha = 1 fit window
frac_latefit   = [0.50 0.70];      % free alpha fit window
min_fit_pts    = 3;                % minimum points required in a fit window
% --------------------------------------------------------------------

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
fprintf('MSD mean type: %s\n', msd_mean_type);
fprintf('Loaded %d trajectories\n', nTraj);
fprintf('Trajectory IDs: %s\n', mat2str(used_ids'));
fprintf('Common resampled dt: %.6g s\n', common_dt);
fprintf('Native dt range: [%.6g, %.6g] s\n', min(native_dts), max(native_dts));

%%%% ENSEMBLE MSD UP TO LONGEST AVAILABLE LAG
maxLag = max(track_lengths) - 1;

tau = (1:maxLag)' * common_dt;
msd_sum    = zeros(maxLag,1);
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

        msd_sum(lag)    = msd_sum(lag) + sum(sqdisp);
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

% only use valid points for fitting
valid = isfinite(ensemble_msd) & ensemble_msd > 0 & pair_count > 0;
tau_valid = tau(valid);
msd_valid = ensemble_msd(valid);
nValid = numel(tau_valid);

if nValid < min_fit_pts
    error('Not enough valid MSD points to fit.');
end

%%%% HELPERS FOR FRACTIONAL WINDOWS
to_idx = @(fr) max(1, min(nValid, round(fr * nValid)));

idx_ball = to_idx(frac_ballistic(1)) + 1 : to_idx(frac_ballistic(2));
idx_diff = to_idx(frac_diffusive(1)) + 1 : to_idx(frac_diffusive(2));
idx_late = to_idx(frac_latefit(1))   + 1 : to_idx(frac_latefit(2));

if numel(idx_ball) < min_fit_pts || numel(idx_diff) < min_fit_pts || numel(idx_late) < min_fit_pts
    error('One or more fit windows contain too few points. Adjust the fraction settings.');
end

%%%% FITS ON LOG-LOG AXES
logtau = log10(tau_valid);
logmsd = log10(msd_valid);

% alpha = 2 fixed slope, fit intercept only
b_ball = mean(logmsd(idx_ball) - 2*logtau(idx_ball));
fit_ball = 10.^(2*logtau(idx_ball) + b_ball);

% alpha = 1 fixed slope, fit intercept only
b_diff = mean(logmsd(idx_diff) - 1*logtau(idx_diff));
fit_diff = 10.^(1*logtau(idx_diff) + b_diff);

% free slope fit on late window
p_late = polyfit(logtau(idx_late), logmsd(idx_late), 1);
alpha_late = p_late(1);
K_late = 10^(p_late(2));
fit_late = 10.^(polyval(p_late, logtau(idx_late)));

fprintf('Late-window fitted alpha = %.3f\n', alpha_late);
fprintf('Late-window fitted K = %.3e m^2 / s^alpha\n', K_late);

%%%% SECOND PANEL: TRAJECTORY COUNT PER LAG
traj_count = zeros(maxLag,1);
for lag = 1:maxLag
    c = 0;
    for k = 1:nTraj
        if track_lengths(k) > lag
            c = c + 1;
        end
    end
    traj_count(lag) = c;
end

%%% PLOT: TWO SUBPLOTS WITH MATCHED X-LIMITS
%figure('Color','w');

ax1 = subplot(2,1,1);
loglog(tau, ensemble_msd, 'o-', 'LineWidth', 2, 'DisplayName',cond_name +", micro" );
hold on;
% loglog(tau_valid(idx_ball), fit_ball, '--', 'LineWidth', 2, ...
%     'DisplayName', sprintf('\\alpha = 2 fit (%.0f-%.0f%%%)', 100*frac_ballistic(1), 100*frac_ballistic(2)));
% loglog(tau_valid(idx_diff), fit_diff, '--', 'LineWidth', 2, ...
%     'DisplayName', sprintf('\\alpha = 1 fit (%.0f-%.0f%%%)', 100*frac_diffusive(1), 100*frac_diffusive(2)));
% loglog(tau_valid(idx_late), fit_late, '--', 'LineWidth', 2, ...
%     'DisplayName', sprintf('\\alpha = %.2f fit (%.0f-%.0f%%%)', alpha_late, 100*frac_latefit(1), 100*frac_latefit(2)));
grid on;
xlabel('\tau (s)');
ylabel('Ensemble MSD(\tau) (m^2)');
%title(sprintf('Ensemble MSD - %s', cond_name));
title('Ensemble MSD Comparison, microalgae')
legend('Location','best');

% ax2 = subplot(2,1,2);
% loglog(tau, pair_count, 'o-', 'LineWidth', 2, 'DisplayName', 'Pair count');
% hold on;
% %loglog(tau, traj_count, 's-', 'LineWidth', 1.5, 'DisplayName', 'Trajectory count');
% grid on;
% xlabel('\tau (s)');
% ylabel('Count');
% title('Pair count per lag');
% %legend('Location','best');

rel_err = nan(maxLag,1);
Meff = nan(maxLag,1);

for lag = 1:maxLag
    idx = track_lengths > lag;
    Meff(lag) = sum(1.5 * track_lengths(idx) / lag);
    rel_err(lag) = 1 / sqrt(Meff(lag));   % for 2D diffusion
end

percent_err = 100 * rel_err;

ax2 = subplot(2,1,2);
grid on; hold on;


plot(tau,percent_err,'-',LineWidth=2.5)
ylabel('Error %')
title('Error quantification')
set(gca,"Xscale",'log')

% Match x-axis only
linkaxes([ax1, ax2], 'x');
xlim(ax1, [min(tau_valid) max(tau_valid)]);

% %%%% OPTIONAL SEPARATE FIGURE FOR COUNTS ONLY
% figure('Color','w');
% loglog(tau, pair_count, 'o-', 'LineWidth', 2, 'DisplayName', 'Pair count');
% hold on;
% loglog(tau, traj_count, 's-', 'LineWidth', 1.5, 'DisplayName', 'Trajectory count');
% grid on;
% xlabel('\tau (s)');
% ylabel('Count');
% title(sprintf('Counts per lag - %s', cond_name));
% legend('Location','best');

%%%% SUMMARY
fprintf('Shortest trajectory: %d frames\n', min(track_lengths));
fprintf('Longest trajectory: %d frames\n', max(track_lengths));
fprintf('Max lag analyzed: %d resampled steps = %.3f s\n', maxLag, maxLag*common_dt);
fprintf('Ballistic fit window: %d points\n', numel(idx_ball));
fprintf('Diffusive fit window: %d points\n', numel(idx_diff));
fprintf('Late fit window: %d points\n', numel(idx_late));

% set(ax1,"XScale","linear")
% set(ax2,"XScale","linear")
% set(ax1,"YScale","linear")
% set(ax2,"YScale","linear")

set(gcf,"Position",[100,100,1400,900])
subplot(2,1,1); set(gca,'FontSize',18)
subplot(2,1,2); set(gca,'FontSize',18)
%subplot(2,1,3); set(gca,'FontSize',15)
