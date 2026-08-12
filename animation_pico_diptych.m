clear; clc; close all;

%%%% SETTINGS
dir_w = "./data_standardized/pico_water/";
dir_m = "./data_standardized/pico_mucus/";

traj_ids = 1:36;

go = 3;                             % sgolay order
gw = 11;                            % sgolay window
use_smoothing = true;

min_timesteps_keep = 1;             % discard short tracks if desired
show_trails = true;                 % keep past path visible
trail_length = inf;                 % use inf for full trail, or e.g. 40 for recent tail only

frame_delay = 0.08;                 % seconds between gif frames
gif_name = 'picoalgae_water_vs_mucus.gif';

warning('off','all');

%%%% LOAD WATER TRAJECTORIES
traj_data_w = {};
traj_time_w = {};
track_lengths_w = [];
used_ids_w = [];

for n = traj_ids
    filename = dir_w + sprintf('%i.csv', n);

    if ~isfile(filename)
        fprintf('Skipping missing water file: %s\n', filename);
        continue;
    end

    try
        [x, y, t, ~, meta] = load_standardized_trajectory(filename);
    catch ME
        fprintf('Skipping unreadable water file: %s (%s)\n', filename, ME.message);
        continue;
    end

    if ~meta.has_time
        fprintf('Skipping water file with missing Time_s column: %s\n', filename);
        continue;
    end

    if numel(x) < min_timesteps_keep
        fprintf('Discarding water trajectory %d (%d timesteps)\n', n, numel(x));
        continue;
    end

    if use_smoothing && numel(x) >= gw
        x = sgolayfilt(x, go, gw);
        y = sgolayfilt(y, go, gw);
    end

    traj_data_w{end+1,1} = [x(:), y(:)];
    traj_time_w{end+1,1} = t(:);
    track_lengths_w(end+1,1) = numel(x);
    used_ids_w(end+1,1) = n;
end

%%%% LOAD MUCUS TRAJECTORIES
traj_data_m = {};
traj_time_m = {};
track_lengths_m = [];
used_ids_m = [];

for n = traj_ids
    filename = dir_m + sprintf('%i.csv', n);

    if ~isfile(filename)
        fprintf('Skipping missing mucus file: %s\n', filename);
        continue;
    end

    try
        [x, y, t, ~, meta] = load_standardized_trajectory(filename);
    catch ME
        fprintf('Skipping unreadable mucus file: %s (%s)\n', filename, ME.message);
        continue;
    end

    if ~meta.has_time
        fprintf('Skipping mucus file with missing Time_s column: %s\n', filename);
        continue;
    end

    if numel(x) < min_timesteps_keep
        fprintf('Discarding mucus trajectory %d (%d timesteps)\n', n, numel(x));
        continue;
    end

    if use_smoothing && numel(x) >= gw
        x = sgolayfilt(x, go, gw);
        y = sgolayfilt(y, go, gw);
    end

    traj_data_m{end+1,1} = [x(:), y(:)];
    traj_time_m{end+1,1} = t(:);
    track_lengths_m(end+1,1) = numel(x);
    used_ids_m(end+1,1) = n;
end

warning('on','all');

nTraj_w = numel(traj_data_w);
nTraj_m = numel(traj_data_m);

if nTraj_w == 0 && nTraj_m == 0
    error('No trajectories available after filtering.');
end

maxFrames_w = 0;
maxFrames_m = 0;

if nTraj_w > 0
    maxFrames_w = max(track_lengths_w);
end
if nTraj_m > 0
    maxFrames_m = max(track_lengths_m);
end

maxFrames = max(maxFrames_w, maxFrames_m);

fprintf('Water trajectories kept: %d\n', nTraj_w);
fprintf('Mucus trajectories kept: %d\n', nTraj_m);
fprintf('Max frames (global): %d\n', maxFrames);

%%%% GLOBAL AXIS LIMITS - WATER
if nTraj_w > 0
    allx_w = [];
    ally_w = [];
    for k = 1:nTraj_w
        allx_w = [allx_w; traj_data_w{k}(:,1)];
        ally_w = [ally_w; traj_data_w{k}(:,2)];
    end

    xmin_w = min(allx_w); xmax_w = max(allx_w);
    ymin_w = min(ally_w); ymax_w = max(ally_w);

    xr_w = xmax_w - xmin_w;
    yr_w = ymax_w - ymin_w;
    pad_w = 0.05 * max([xr_w, yr_w, eps]);

    xmin_w = xmin_w - pad_w; xmax_w = xmax_w + pad_w;
    ymin_w = ymin_w - pad_w; ymax_w = ymax_w + pad_w;
else
    xmin_w = -1; xmax_w = 1;
    ymin_w = -1; ymax_w = 1;
end

%%%% GLOBAL AXIS LIMITS - MUCUS
if nTraj_m > 0
    allx_m = [];
    ally_m = [];
    for k = 1:nTraj_m
        allx_m = [allx_m; traj_data_m{k}(:,1)];
        ally_m = [ally_m; traj_data_m{k}(:,2)];
    end

    xmin_m = min(allx_m); xmax_m = max(allx_m);
    ymin_m = min(ally_m); ymax_m = max(ally_m);

    xr_m = xmax_m - xmin_m;
    yr_m = ymax_m - ymin_m;
    pad_m = 0.05 * max([xr_m, yr_m, eps]);

    xmin_m = xmin_m - pad_m; xmax_m = xmax_m + pad_m;
    ymin_m = ymin_m - pad_m; ymax_m = ymax_m + pad_m;
else
    xmin_m = -1; xmax_m = 1;
    ymin_m = -1; ymax_m = 1;
end

%%%% COLORS
colors_w = lines(max(nTraj_w,1));
colors_m = lines(max(nTraj_m,1));

%%%% FIGURE SETUP
fig = figure('Color','w','Position',[100 100 1400 700]);

ax1 = subplot(1,2,1, 'Parent', fig);
hold(ax1, 'on');
axis(ax1, 'equal');
xlim(ax1, [xmin_w xmax_w]);
ylim(ax1, [ymin_w ymax_w]);
grid(ax1, 'on');
xlabel(ax1, 'x (m)');
ylabel(ax1, 'y (m)');
title(ax1, 'Picoalgae trajectories in water');

ax2 = subplot(1,2,2, 'Parent', fig);
hold(ax2, 'on');
axis(ax2, 'equal');
xlim(ax2, [xmin_m xmax_m]);
ylim(ax2, [ymin_m ymax_m]);
grid(ax2, 'on');
xlabel(ax2, 'x (m)');
ylabel(ax2, 'y (m)');
title(ax2, 'Picoalgae trajectories in mucus');
axis equal

%%%% PRECREATE GRAPHICS OBJECTS - WATER
trail_handles_w = gobjects(nTraj_w,1);
head_handles_w  = gobjects(nTraj_w,1);

for k = 1:nTraj_w
    trail_handles_w(k) = plot(ax1, nan, nan, '-', 'LineWidth', 1.5, ...
        'Color', colors_w(k,:));
    head_handles_w(k) = plot(ax1, nan, nan, 'o', ...
        'MarkerSize', 6, ...
        'MarkerFaceColor', colors_w(k,:), ...
        'MarkerEdgeColor', colors_w(k,:));
end

%%%% PRECREATE GRAPHICS OBJECTS - MUCUS
trail_handles_m = gobjects(nTraj_m,1);
head_handles_m  = gobjects(nTraj_m,1);

for k = 1:nTraj_m
    trail_handles_m(k) = plot(ax2, nan, nan, '-', 'LineWidth', 1.5, ...
        'Color', colors_m(k,:));
    head_handles_m(k) = plot(ax2, nan, nan, 'o', ...
        'MarkerSize', 6, ...
        'MarkerFaceColor', colors_m(k,:), ...
        'MarkerEdgeColor', colors_m(k,:));
end

time_text_w = text(ax1, 0.02, 0.98, '', 'Units', 'normalized', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 12, 'FontWeight', 'bold');

time_text_m = text(ax2, 0.02, 0.98, '', 'Units', 'normalized', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 12, 'FontWeight', 'bold');

%%%% ANIMATE AND WRITE GIF
for f = 1:maxFrames

    % --- WATER ---
    for k = 1:nTraj_w
        xy = traj_data_w{k};
        N = size(xy,1);

        if f <= N
            if show_trails
                if isinf(trail_length)
                    idx1 = 1;
                else
                    idx1 = max(1, f - trail_length + 1);
                end
                idx2 = f;
                set(trail_handles_w(k), 'XData', xy(idx1:idx2,1), 'YData', xy(idx1:idx2,2));
            else
                set(trail_handles_w(k), 'XData', nan, 'YData', nan);
            end

            set(head_handles_w(k), 'XData', xy(f,1), 'YData', xy(f,2));
        else
            if ~show_trails
                set(trail_handles_w(k), 'XData', nan, 'YData', nan);
            end
            set(head_handles_w(k), 'XData', nan, 'YData', nan);
        end
    end

    % --- MUCUS ---
    for k = 1:nTraj_m
        xy = traj_data_m{k};
        N = size(xy,1);

        if f <= N
            if show_trails
                if isinf(trail_length)
                    idx1 = 1;
                else
                    idx1 = max(1, f - trail_length + 1);
                end
                idx2 = f;
                set(trail_handles_m(k), 'XData', xy(idx1:idx2,1), 'YData', xy(idx1:idx2,2));
            else
                set(trail_handles_m(k), 'XData', nan, 'YData', nan);
            end

            set(head_handles_m(k), 'XData', xy(f,1), 'YData', xy(f,2));
        else
            if ~show_trails
                set(trail_handles_m(k), 'XData', nan, 'YData', nan);
            end
            set(head_handles_m(k), 'XData', nan, 'YData', nan);
        end
    end

    if nTraj_w > 0
        current_time_w = max(cellfun(@(t) t(min(f, numel(t))), traj_time_w));
        total_time_w = max(cellfun(@(t) t(end), traj_time_w));
        set(time_text_w, 'String', sprintf('Time: %.2f / %.2f s', current_time_w, total_time_w));
    end

    if nTraj_m > 0
        current_time_m = max(cellfun(@(t) t(min(f, numel(t))), traj_time_m));
        total_time_m = max(cellfun(@(t) t(end), traj_time_m));
        set(time_text_m, 'String', sprintf('Time: %.2f / %.2f s', current_time_m, total_time_m));
    end

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
