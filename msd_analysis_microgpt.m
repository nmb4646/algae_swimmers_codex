clear; clc; close all;

blue = [.06,.43,.74];

dir_w = './data/algae_motion/buffer_micro/';
dir_m = './data/algae_motion/mucus_micro/';

% ---------------- USER SETTINGS ----------------
use_mucus = true;          % true = mucus, false = water
traj_ids = 1:20;

dt  = 1/16.5;             % s
pix = 0.24e-6;            % m / pixel

go = 3;                   % sgolay order
gw = 11;                  % sgolay window

min_timesteps_keep = 70; % discard trajectories shorter than this
min_pairs = 20;           % minimum displacement pairs for reliable MSD point
alpha_win = 3;            % sliding window size for local alpha (odd recommended)
late_frac = 0.5;          % fraction of valid lag points used for "late-time" alpha
min_late_pts = 50;         % at least this many points in late-time fit
% ------------------------------------------------

if use_mucus
    dir_data = dir_m;
    file_pattern = "Microalgae in mucus (%i).csv";
    cond_name = "mucus";
else
    dir_data = dir_w;
    file_pattern = "Microalgae in water (%i).csv";
    cond_name = "water";
end

% Storage
traj_data = {};
track_lengths = [];
used_ids = [];
discarded_ids = [];

%%% LOAD AND FILTER TRAJECTORIES
warning('off','all');

for n = traj_ids
    filename = dir_data + sprintf(file_pattern, n);

    if ~isfile(filename)
        fprintf('Skipping missing file: %s\n', filename);
        discarded_ids(end+1) = n;
        continue;
    end

    T = readtable(filename);

    if ~all(ismember({'X','Y'}, T.Properties.VariableNames))
        fprintf('Skipping file with missing X/Y columns: %s\n', filename);
        discarded_ids(end+1) = n;
        continue;
    end

    x = T.X;
    y = T.Y;
    
    %[x,y] = repeat_traj_rot4(x,y);
    maxT = length(x);

    

    if maxT < min_timesteps_keep
        fprintf('Discarding trajectory %d (%d timesteps < %d)\n', n, maxT, min_timesteps_keep);
        discarded_ids(end+1) = n;
        continue;
    end

    % shift to origin and convert to meters
    x = x - x(1);
    y = y - y(1);
    x = x * pix;
    y = y * pix;

    % smooth
    x = sgolayfilt(x, go, gw);
    y = sgolayfilt(y, go, gw);

    traj_data{end+1,1} = [x(:), y(:)];
    track_lengths(end+1,1) = numel(x);
    used_ids(end+1,1) = n;
end

warning('on','all');

nTraj = numel(traj_data);

if nTraj == 0
    error('No trajectories passed the min_timesteps_keep threshold.');
end

fprintf('\nCondition: %s\n', cond_name);
fprintf('Trajectories kept: %d\n', nTraj);
fprintf('Trajectory IDs kept: %s\n', mat2str(used_ids'));
fprintf('Trajectory IDs discarded: %s\n', mat2str(discarded_ids));
fprintf('Shortest kept trajectory: %d timesteps\n', min(track_lengths));
fprintf('Longest kept trajectory: %d timesteps\n\n', max(track_lengths));

%%% ENSEMBLE MSD
% Max lag limited only by shortest KEPT trajectory
maxLag = min(track_lengths) - 1;
tau = (1:maxLag)' * dt;

msd_sum   = zeros(maxLag,1);  % sum of all squared displacements
msd_count = zeros(maxLag,1);  % number of displacement pairs
msd_each  = nan(nTraj, maxLag);

for k = 1:nTraj
    xy = traj_data{k};
    x = xy(:,1);
    y = xy(:,2);
    N = numel(x);

    for lag = 1:maxLag
        dx = x(1+lag:N) - x(1:N-lag);
        dy = y(1+lag:N) - y(1:N-lag);
        sqdisp = dx.^2 + dy.^2;

        msd_each(k,lag) = mean(sqdisp);
        msd_sum(lag)    = msd_sum(lag) + sum(sqdisp);
        msd_count(lag)  = msd_count(lag) + numel(sqdisp);
    end
end

ensemble_msd = msd_sum ./ msd_count;
traj_count = sum(~isnan(msd_each), 1)';

%%% LOCAL ANOMALOUS EXPONENT alpha(tau)
% Keep only lag points with enough displacement pairs and positive MSD
valid = (msd_count >= min_pairs) & isfinite(ensemble_msd) & (ensemble_msd > 0);

tau_valid = tau(valid);
msd_valid = ensemble_msd(valid);

if numel(tau_valid) < 3
    error('Too few valid lag points after applying min_pairs threshold.');
end

logtau = log10(tau_valid);
logmsd = log10(msd_valid);

halfwin = floor(alpha_win/2);
alpha_tau = nan(size(tau_valid));
alpha_ci = nan(numel(tau_valid),2);

for i = 1:numel(tau_valid)
    i1 = max(1, i-halfwin);
    i2 = min(numel(tau_valid), i+halfwin);

    if (i2 - i1 + 1) >= 3
        p = polyfit(logtau(i1:i2), logmsd(i1:i2), 1);
        alpha_tau(i) = p(1);

        mdl = fitlm(logtau(i1:i2), logmsd(i1:i2));
        ci = coefCI(mdl, 0.05);  % intercept row 1, slope row 2
        alpha_ci(i,:) = ci(2,:);
    end
end

%%% LATE-TIME ALPHA
nValid = numel(tau_valid);
nLate = max(min_late_pts, round(late_frac * nValid));
nLate = min(nLate, nValid);
idxLate = (nValid - nLate + 1):nValid;

pLate = polyfit(logtau(idxLate), logmsd(idxLate), 1);
alpha_late = pLate(1);
Kalpha_late = 10^(pLate(2));   % MSD = K_alpha * tau^alpha

mdlLate = fitlm(logtau(idxLate), logmsd(idxLate));
ciLate = coefCI(mdlLate, 0.05);
alpha_late_ci = ciLate(2,:);

%%% OPTIONAL: SIMPLE PER-TRAJECTORY METRICS FROM ORIGINAL STYLE
v_mean = nan(nTraj,1);
omega_mean = nan(nTraj,1);

for k = 1:nTraj
    xy = traj_data{k};
    x = xy(:,1);
    y = xy(:,2);

    dx = gradient(x)/dt;
    dy = gradient(y)/dt;
    speed = sqrt(dx.^2 + dy.^2);

    theta = atan2(dy, dx);
    theta_unwrapped = unwrap(theta);
    total_turn = theta_unwrapped(end) - theta_unwrapped(1);

    omega_mean(k) = total_turn/(length(x)*dt);
    v_mean(k) = mean(speed);
end

%%% PLOTS
figure;
loglog(tau, ensemble_msd, 'o-', 'LineWidth', 2, 'Color', blue, ...
    'DisplayName', 'Ensemble MSD');
hold on;
loglog(tau_valid(idxLate), 10.^polyval(pLate, log10(tau_valid(idxLate))), '--k', ...
    'LineWidth', 2, ...
    'DisplayName', sprintf('Late fit: \\alpha = %.2f', alpha_late));
xlabel('\tau (s)');
ylabel('MSD(\tau) (m^2)');
title(sprintf('Ensemble MSD (%s)', cond_name));
legend('Location','best');
grid on;

figure;
semilogy(tau, msd_count, 'o-', 'LineWidth', 2, 'Color', [0.85 0.33 0.10]);
xlabel('\tau (s)');
ylabel('Number of displacement pairs');
title(sprintf('MSD pair counts (%s)', cond_name));
grid on;

figure;
plot(tau, traj_count, 'o-', 'LineWidth', 2, 'Color', [0.2 0.6 0.2]);
xlabel('\tau (s)');
ylabel('Number of trajectories contributing');
title(sprintf('Trajectory counts per lag (%s)', cond_name));
grid on;

figure;
plot(tau_valid, alpha_tau, 'o-', 'LineWidth', 2, 'Color', [0.49 0.18 0.56]);
hold on;
yline(2, '--', 'Ballistic');
yline(1, '--', 'Diffusive');
yline(0, '--', 'Confined');
xlabel('\tau (s)');
ylabel('\alpha(\tau)');
title(sprintf('Local anomalous exponent (%s)', cond_name));
grid on;

%%% DISPLAY SUMMARY
fprintf('========== SUMMARY ==========\n');
fprintf('Condition: %s\n', cond_name);
fprintf('Kept %d trajectories\n', nTraj);
fprintf('maxLag = %d frames (%.3f s)\n', maxLag, maxLag*dt);
fprintf('Valid lag points after min_pairs threshold: %d / %d\n', numel(tau_valid), numel(tau));
fprintf('Late-time alpha = %.3f\n', alpha_late);
fprintf('95%%% CI for late-time alpha = [%.3f, %.3f]\n', alpha_late_ci(1), alpha_late_ci(2));
fprintf('Late-time K_alpha = %.3e [m^2 / s^alpha]\n', Kalpha_late);
fprintf('Mean speed across kept trajectories = %.3e m/s\n', mean(v_mean,'omitnan'));
fprintf('Mean angular velocity across kept trajectories = %.3e rad/s\n', mean(omega_mean,'omitnan'));
fprintf('=============================\n');