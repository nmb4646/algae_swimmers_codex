clear; clc; close all;

%%% SETTINGS
dir_m = './data/algae_motion/mucus_micro/';
n = 9;                          % trajectory number to analyze

dt  = 1/16.5;                   % s
pix = 0.24e-6;                  % m / pixel

go = 3;                         % sgolay order
gw = 11;                        % sgolay window

% ---------- trap detection parameters ----------
speed_win_sec = 1.5;              % rolling window for speed averaging
conf_win_sec  = 2.0;            % rolling window for confinement check
min_trap_sec  = 2.0;            % must remain trapped for at least this long

speed_thresh  = 6e-6;          % m/s ; adjust for your data
radius_thresh = 3e-6;           % m   ; max local radius to count as trapped
% ----------------------------------------------

%%% LOAD TRAJECTORY
filename = dir_m + sprintf("Microalgae in mucus (%i).csv", n);

warning('off','all');
T = readtable(filename);
warning('on','all');

x = (T.X - T.X(1)) * pix;
y = (T.Y - T.Y(1)) * pix;

if numel(x) < gw
    error('Trajectory is shorter than Savitzky-Golay window.');
end

% Smooth
x = sgolayfilt(x, go, gw);
y = sgolayfilt(y, go, gw);

N = numel(x);
t = (0:N-1)' * dt;

%%% INSTANTANEOUS SPEED
dx = gradient(x) / dt;
dy = gradient(y) / dt;
speed = sqrt(dx.^2 + dy.^2);

%speed_thresh=.5*max(speed);

%%% ROLLING SPEED
speed_win = max(3, round(speed_win_sec / dt));
if mod(speed_win,2)==0
    speed_win = speed_win + 1;
end

speed_roll = movmean(speed, speed_win, 'Endpoints','shrink');

%%% ROLLING CONFINEMENT RADIUS
% For each frame, compute the radius of points in a local window
% around their local centroid. Small radius => confined/caged.
conf_win = max(3, round(conf_win_sec / dt));
if mod(conf_win,2)==0
    conf_win = conf_win + 1;
end

half_conf = floor(conf_win/2);
local_radius = nan(N,1);

for i = 1:N
    i1 = max(1, i-half_conf);
    i2 = min(N, i+half_conf);

    xw = x(i1:i2);
    yw = y(i1:i2);

    xc = mean(xw);
    yc = mean(yw);

    r = sqrt((xw-xc).^2 + (yw-yc).^2);
    local_radius(i) = max(r);   % could also use mean(r) or rms(r)
end

%%% TRAP LOGIC
slow_mask   = speed_roll < speed_thresh;
conf_mask   = local_radius < radius_thresh;
trap_mask   = slow_mask & conf_mask;

% Require the mask to stay true for a minimum duration
min_trap_frames = max(3, round(min_trap_sec / dt));

trap_start_idx = NaN;
trap_end_idx   = NaN;

runlen = 0;
for i = 1:N
    if trap_mask(i)
        runlen = runlen + 1;
    else
        runlen = 0;
    end

    if runlen >= min_trap_frames
        trap_start_idx = i - min_trap_frames + 1;

        % extend to end of that continuous trapped run
        j = i;
        while j < N && trap_mask(j+1)
            j = j + 1;
        end
        trap_end_idx = j;
        break
    end
end

%%% RESULTS
if isnan(trap_start_idx)
    fprintf('No trapped regime detected for trajectory %d.\n', n);
else
    fprintf('Trajectory %d trapped starting at frame %d (t = %.2f s)\n', ...
        n, trap_start_idx, t(trap_start_idx));
    fprintf('Trapped segment ends at frame %d (t = %.2f s)\n', ...
        trap_end_idx, t(trap_end_idx));
end

%%% PLOTS
figure('Color','w','Position',[100 100 1200 900]);

% Trajectory
subplot(3,1,1); hold on;
plot(x*1e6, y*1e6, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
plot(x(1)*1e6, y(1)*1e6, 'go', 'MarkerFaceColor', 'g');
plot(x(end)*1e6, y(end)*1e6, 'ro', 'MarkerFaceColor', 'r');

if ~isnan(trap_start_idx)
    plot(x(trap_start_idx:end)*1e6, y(trap_start_idx:end)*1e6, ...
        'b-', 'LineWidth', 2.5);
    plot(x(trap_start_idx)*1e6, y(trap_start_idx)*1e6, ...
        'bo', 'MarkerFaceColor', 'b');
end

axis equal;
grid on;
xlabel('x (\mum)');
ylabel('y (\mum)');
title(sprintf('Trajectory %d', n));
legend('Full path','Start','End','Detected trapped regime','Trap onset','Location','best');

% Speed
subplot(3,1,2); hold on;
plot(t, speed*1e6, '-', 'Color', [0.8 0.8 0.8]);
plot(t, speed_roll*1e6, 'k-', 'LineWidth', 2);
yline(speed_thresh*1e6, '--r', 'Speed threshold');

if ~isnan(trap_start_idx)
    xline(t(trap_start_idx), '-b', 'Trap onset', 'LineWidth', 2);
    xline(t(trap_end_idx), '--b', 'Trap end', 'LineWidth', 1.5);
end

grid on;
xlabel('Time (s)');
ylabel('Speed (\mum/s)');
title('Instantaneous and rolling speed');

% Confinement radius
subplot(3,1,3); hold on;
plot(t, local_radius*1e6, 'm-', 'LineWidth', 2);
yline(radius_thresh*1e6, '--r', 'Radius threshold');

if ~isnan(trap_start_idx)
    xline(t(trap_start_idx), '-b', 'Trap onset', 'LineWidth', 2);
    xline(t(trap_end_idx), '--b', 'Trap end', 'LineWidth', 1.5);
end

grid on;
xlabel('Time (s)');
ylabel('Local radius (\mum)');
title('Local confinement metric');

%%% OPTIONAL: return trapped segment only
if ~isnan(trap_start_idx)
    x_trapped = x(trap_start_idx:trap_end_idx);
    y_trapped = y(trap_start_idx:trap_end_idx);
else
    x_trapped = [];
    y_trapped = [];
end