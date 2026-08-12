clear; clc; close all;

%%% SETTINGS
dir_w = './data/algae_motion/buffer/';
dir_m = './data/algae_motion/mucus/';

use_mucus = false;                    % true = mucus, false = water
traj_ids = 1:36;

pix = 0.24e-6;                      % m / pixel
go = 3;                             % sgolay order
gw = 11;                            % sgolay window
use_smoothing = true;

min_timesteps_keep = 1;             % discard short tracks if desired
show_trails = true;                 % keep past path visible
trail_length = inf;                 % use inf for full trail, or e.g. 40 for recent tail only

frame_delay = 0.08;                 % seconds between gif frames
gif_name = 'microalgae_all_tracks.gif';

%%% FILE PATTERN
if use_mucus
    dir_data = dir_m;
    file_pattern = "%i.csv";
    cond_name = "mucus";
else
    dir_data = dir_w;
    file_pattern = "%i.csv";
    cond_name = "water";
end

%%% LOAD TRAJECTORIES
traj_data = {};
track_lengths = [];
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
        fprintf('Skipping file with missing X/Y columns: %s\n', filename);
        continue;
    end

    x = (T.X - T.X(1)) * pix;
    y = (T.Y - T.Y(1)) * pix;

    if numel(x) < min_timesteps_keep
        fprintf('Discarding trajectory %d (%d timesteps)\n', n, numel(x));
        continue;
    end

    if use_smoothing && numel(x) >= gw
        x = sgolayfilt(x, go, gw);
        y = sgolayfilt(y, go, gw);
    end

    traj_data{end+1,1} = [x(:), y(:)];
    track_lengths(end+1,1) = numel(x);
    used_ids(end+1,1) = n;
end

warning('on','all');

nTraj = numel(traj_data);

if nTraj == 0
    error('No trajectories available after filtering.');
end

maxFrames = max(track_lengths);

fprintf('Condition: %s\n', cond_name);
fprintf('Trajectories kept: %d\n', nTraj);
fprintf('Max frames: %d\n', maxFrames);

%%% GLOBAL AXIS LIMITS
allx = [];
ally = [];
for k = 1:nTraj
    allx = [allx; traj_data{k}(:,1)];
    ally = [ally; traj_data{k}(:,2)];
end

xmin = min(allx); xmax = max(allx);
ymin = min(ally); ymax = max(ally);

% Add padding
xr = xmax - xmin;
yr = ymax - ymin;
pad = 0.05 * max([xr, yr, eps]);

xmin = xmin - pad; xmax = xmax + pad;
ymin = ymin - pad; ymax = ymax + pad;

%%% COLORS
colors = lines(nTraj);

%%% FIGURE SETUP
fig = figure('Color','w','Position',[100 100 800 800]);
ax = axes(fig);
hold(ax, 'on');
axis(ax, 'equal');
xlim(ax, [xmin xmax]);
ylim(ax, [ymin ymax]);
grid(ax, 'on');

xlabel(ax, 'x (m)');
ylabel(ax, 'y (m)');
title(ax, sprintf('Microalgae trajectories in %s', cond_name));

% Precreate graphics objects
trail_handles = gobjects(nTraj,1);
head_handles  = gobjects(nTraj,1);

for k = 1:nTraj
    trail_handles(k) = plot(ax, nan, nan, '-', 'LineWidth', 1.5, ...
        'Color', colors(k,:));
    head_handles(k) = plot(ax, nan, nan, 'o', ...
        'MarkerSize', 6, ...
        'MarkerFaceColor', colors(k,:), ...
        'MarkerEdgeColor', colors(k,:));
end

time_text = text(ax, 0.02, 0.98, '', 'Units', 'normalized', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 12, 'FontWeight', 'bold');

%%% ANIMATE AND WRITE GIF
for f = 1:maxFrames
    for k = 1:nTraj
        xy = traj_data{k};
        N = size(xy,1);

        if f <= N
            if show_trails
                if isinf(trail_length)
                    idx1 = 1;
                else
                    idx1 = max(1, f - trail_length + 1);
                end
                idx2 = f;
                set(trail_handles(k), 'XData', xy(idx1:idx2,1), 'YData', xy(idx1:idx2,2));
            else
                set(trail_handles(k), 'XData', nan, 'YData', nan);
            end

            set(head_handles(k), 'XData', xy(f,1), 'YData', xy(f,2));
        else
            % Keep final trail, hide head after track ends
            if ~show_trails
                set(trail_handles(k), 'XData', nan, 'YData', nan);
            end
            set(head_handles(k), 'XData', nan, 'YData', nan);
        end
    end

    set(time_text, 'String', sprintf('Frame: %d / %d', f, maxFrames));

    drawnow;

    frame = getframe(fig);
    im = frame2im(frame);
    [A, map] = rgb2ind(im, 256);

    if f == 1
        imwrite(A, map, gif_name, 'gif', 'LoopCount', inf, 'DelayTime', frame_delay);
    else
        imwrite(A, map, gif_name, 'gif', 'WriteMode', 'append', 'DelayTime', frame_delay);
    end
end

fprintf('GIF saved as: %s\n', gif_name);