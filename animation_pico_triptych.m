clear; clc; close all;

%%%% SETTINGS
base_dir = "./data/Picoalgae in 0-2% mucin/";

condition_dirs = [ ...
    "picoalgae 0% mucin"; ...
    "picoalgae 1% mucin"; ...
    "picoalgae 2% mucin"];

condition_names = ["0% mucin", "1% mucin", "2% mucin"];

fps = 15;
dt = 1 / fps;
pix_m_per_px = 0.48e-6;

go = 3;                             % sgolay order
gw = 11;                            % sgolay window
use_smoothing = true;

min_timesteps_keep = 100;             % discard short tracks if desired
show_trails = true;                 % keep past path visible
trail_length = inf;                 % use inf for full trail, or e.g. 40 for recent tail only

frame_delay = 0.08;                 % seconds between gif frames
gif_name = 'picoalgae_0_1_2_mucin.gif';

if isfile(gif_name)
    delete(gif_name);
end

warning('off','all');

%%%% LOAD TRAJECTORIES
nCond = numel(condition_names);
traj_data = cell(nCond, 1);
traj_time = cell(nCond, 1);
track_lengths = cell(nCond, 1);
used_files = cell(nCond, 1);

for c = 1:nCond
    dir_data = fullfile(base_dir, condition_dirs(c));

    [traj_data{c}, traj_time{c}, track_lengths{c}, used_files{c}] = ...
        load_raw_pico_mucin_animation_trajectories( ...
            dir_data, pix_m_per_px, dt, use_smoothing, go, gw, min_timesteps_keep);
end

warning('on','all');

nTraj = cellfun(@numel, traj_data);

if all(nTraj == 0)
    error('No trajectories available after filtering.');
end

maxFrames = max(cellfun(@max_or_zero, track_lengths));

for c = 1:nCond
    fprintf('%s trajectories kept: %d\n', condition_names(c), nTraj(c));
    fprintf('%s files: %s\n', condition_names(c), strjoin(used_files{c}, ", "));
end
fprintf('Max frames (global): %d\n', maxFrames);

%%%% SHARED AXIS LIMITS
allx_global = [];
ally_global = [];

for c = 1:nCond
    if nTraj(c) > 0
        for k = 1:nTraj(c)
            allx_global = [allx_global; traj_data{c}{k}(:,1)]; %#ok<AGROW>
            ally_global = [ally_global; traj_data{c}{k}(:,2)]; %#ok<AGROW>
        end
    end
end

if isempty(allx_global)
    axis_limits = [-1, 1, -1, 1];
else
    xmin = min(allx_global); xmax = max(allx_global);
    ymin = min(ally_global); ymax = max(ally_global);

    xmid = 0.5 * (xmin + xmax);
    ymid = 0.5 * (ymin + ymax);
    span = max([xmax - xmin, ymax - ymin, eps]);
    pad = 0.05 * span;
    half_span = 0.5 * span + pad;

    axis_limits = [xmid - half_span, xmid + half_span, ...
                   ymid - half_span, ymid + half_span];
end

%%%% COLORS
colors = cell(nCond, 1);
for c = 1:nCond
    colors{c} = lines(max(nTraj(c), 1));
end

%%%% FIGURE SETUP
fig = figure('Color','w','Position',[100 100 1800 650]);

axes_handles = gobjects(nCond, 1);
trail_handles = cell(nCond, 1);
head_handles = cell(nCond, 1);
time_text = gobjects(nCond, 1);

for c = 1:nCond
    axes_handles(c) = subplot(1, 3, c, 'Parent', fig);
    ax = axes_handles(c);
    hold(ax, 'on');
    axis(ax, 'equal');
    pbaspect(ax, [1 1 1]);
    xlim(ax, axis_limits(1:2));
    ylim(ax, axis_limits(3:4));
    % xlim([0 2.5e-4])
    % ylim([0 2.5e-4])
    grid(ax, 'on');
    xlabel(ax, 'x (m)');
    ylabel(ax, 'y (m)');
    title(ax, sprintf('Picoalgae trajectories in %s', condition_names(c)));

    trail_handles{c} = gobjects(nTraj(c), 1);
    head_handles{c} = gobjects(nTraj(c), 1);

    for k = 1:nTraj(c)
        trail_handles{c}(k) = plot(ax, nan, nan, '-', ...
            'LineWidth', 1.5, ...
            'Color', colors{c}(k,:));
        head_handles{c}(k) = plot(ax, nan, nan, 'o', ...
            'MarkerSize', 6, ...
            'MarkerFaceColor', colors{c}(k,:), ...
            'MarkerEdgeColor', colors{c}(k,:));
    end

    time_text(c) = text(ax, 0.02, 0.98, '', 'Units', 'normalized', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', ...
        'FontSize', 12, 'FontWeight', 'bold');
end

%%%% ANIMATE AND WRITE GIF
for f = 1:maxFrames
    for c = 1:nCond
        for k = 1:nTraj(c)
            xy = traj_data{c}{k};
            N = size(xy, 1);

            if f <= N
                if show_trails
                    if isinf(trail_length)
                        idx1 = 1;
                    else
                        idx1 = max(1, f - trail_length + 1);
                    end
                    idx2 = f;
                    set(trail_handles{c}(k), ...
                        'XData', xy(idx1:idx2,1), ...
                        'YData', xy(idx1:idx2,2));
                else
                    set(trail_handles{c}(k), 'XData', nan, 'YData', nan);
                end

                set(head_handles{c}(k), 'XData', xy(f,1), 'YData', xy(f,2));
            else
                if ~show_trails
                    set(trail_handles{c}(k), 'XData', nan, 'YData', nan);
                end
                set(head_handles{c}(k), 'XData', nan, 'YData', nan);
            end
        end

        if nTraj(c) > 0
            current_time = max(cellfun(@(t) t(min(f, numel(t))), traj_time{c}));
            total_time = max(cellfun(@(t) t(end), traj_time{c}));
            set(time_text(c), 'String', sprintf('Time: %.2f / %.2f s', current_time, total_time));
        end
    end

    drawnow;

    if ~isgraphics(fig, 'figure')
        error('Animation figure was closed before the GIF finished writing.');
    end

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

function out = max_or_zero(v)
    if isempty(v)
        out = 0;
    else
        out = max(v);
    end
end

function [traj_data, traj_time, track_lengths, used_files] = ...
    load_raw_pico_mucin_animation_trajectories( ...
        dir_data, pix_m_per_px, dt, use_smoothing, go, gw, min_timesteps_keep)

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
    traj_time = {};
    track_lengths = [];
    used_files = strings(0, 1);

    for k = 1:numel(files)
        filename = fullfile(files(k).folder, files(k).name);

        if dir(filename).bytes == 0
            fprintf('Skipping empty file: %s\n', filename);
            continue;
        end

        try
            T = readtable(filename);
        catch ME
            fprintf('Skipping unreadable file: %s (%s)\n', filename, ME.message);
            continue;
        end

        vars = string(T.Properties.VariableNames);
        if ~all(ismember(["X", "Y"], vars))
            fprintf('Skipping file with missing X/Y columns: %s\n', filename);
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

        if numel(x_px) < min_timesteps_keep
            fprintf('Discarding short trajectory: %s (%d timesteps)\n', filename, numel(x_px));
            continue;
        end

        [frame_number, unique_idx] = unique(frame_number, 'stable');
        x_px = x_px(unique_idx);
        y_px = y_px(unique_idx);

        x = x_px * pix_m_per_px;
        y = y_px * pix_m_per_px;
        t = (frame_number - frame_number(1)) * dt;

        if use_smoothing && numel(x) >= gw
            x = sgolayfilt(x, go, gw);
            y = sgolayfilt(y, go, gw);
        end

        if any(~isfinite(x)) || any(~isfinite(y)) || any(~isfinite(t))
            fprintf('Skipping file with non-finite values after preprocessing: %s\n', filename);
            continue;
        end

        traj_data{end+1,1} = [x(:), y(:)]; %#ok<AGROW>
        traj_time{end+1,1} = t(:); %#ok<AGROW>
        track_lengths(end+1,1) = numel(x); %#ok<AGROW>
        used_files(end+1,1) = string(files(k).name); %#ok<AGROW>
    end
end
