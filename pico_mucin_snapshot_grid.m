clear; clc; close all;

%%%% SETTINGS
% "triptych" reproduces the data source used by animation_pico_triptych.m.
% Change to "latest" to use the updated 8_12_26 delivery instead.
data_source = "triptych";

% "rows"    -> conditions are rows and times are columns (2-by-3)
% "columns" -> conditions are columns and times are rows (3-by-2)
condition_layout = "rows";

snapshot_times_s = [1, 5, 25];
output_png = 'picoalgae_0_2_mucin_snapshots_5_15_25s.png';
output_resolution_dpi = 300;

condition_dirs = [ ...
    "picoalgae 0% mucin"; ...
    "picoalgae 1% mucin"; ...
    "picoalgae 2% mucin"];
condition_names = ["0% mucin", "1% mucin", "2% mucin"];
selected_conditions = [1, 3];       % plot 0% and 2% mucin

% Curated exclusions for the original triptych delivery. Repeat groups are
% listed separately from manually rejected outliers for traceability.
triptych_excluded_files = { ...
    ["7.csv"; ...                    % long dark-red 0% trajectory
     "11.csv"]; ...                  % anomalously straight 0% trajectory
    strings(0,1); ...
    strings(0,1)};
remove_repeat_trajectories = true;
triptych_repeat_files = { ...
    ["5.csv"; ...                    % near-repeat of 3.csv
     "9.csv"; "10.csv"; ...         % repeats/continuation of 6.csv
     "15.csv"; "44.csv"; ...        % truncated repeats of 18.csv
     "16.csv"; ...                   % near-repeat of 20.csv
     "27.csv"; ...                   % three-frame-offset repeat of 28.csv
     "37.csv"; ...                   % frame-offset repeat of 34.csv
     "41.csv"; ...                   % near-repeat of 23.csv
     "42.csv"]; ...                  % truncated repeat of 25.csv
    strings(0,1); ...
    ["14.csv"; ...                   % near-repeat of 1.csv
     "19.csv"; ...                   % near-repeat of 6.csv with a long upward tail
     "23.csv"; ...                   % truncated repeat of 31.csv
     "30.csv"]};                     % near-repeat of 34.csv

fps = 15;
dt = 1 / fps;
pix_mm_per_px = 0.00048;            % documented 0.48 um/px calibration, in mm/px
major_tick_step_mm = 0.25;          % puts grid lines exactly on every box edge
crop_limits_mm = [0, 1, 0, 1];      % [xmin xmax ymin ymax]

% When enabled, translate whole trajectories to reduce intersections and
% near-overlap in the final requested snapshot. Shapes and timing are unchanged.
minimize_trajectory_overlap = true;
overlap_clearance_mm = 0.02;        % desired centerline separation (20 um)
overlap_grid_step_mm = 0.005;       % resolution of the overlap analysis
overlap_candidate_count = 17;       % candidate translations along each axis
overlap_optimization_passes = 3;

go = 3;                             % Savitzky-Golay polynomial order
gw = 11;                            % Savitzky-Golay window
use_smoothing = true;
min_timesteps_keep = 100;
show_trails = true;
trail_length = inf;                 % inf exactly matches the triptych GIF

switch lower(data_source)
    case "triptych"
        base_dir = "./data/Picoalgae in 0-2% mucin/";
        excluded_files = triptych_excluded_files;
        repeat_files = triptych_repeat_files;
    case "latest"
        base_dir = "./data/Picoalgae in 0-2% mucin 8_12_26/New modeling-081026/";
        excluded_files = repmat({strings(0,1)}, numel(condition_names), 1);
        repeat_files = repmat({strings(0,1)}, numel(condition_names), 1);
    otherwise
        error('Unsupported data_source: %s. Use "triptych" or "latest".', data_source);
end

if any(~isfinite(snapshot_times_s)) || any(snapshot_times_s < 0)
    error('snapshot_times_s must contain finite, nonnegative times.');
end

warning('off', 'all');
warning_cleanup = onCleanup(@() warning('on', 'all'));

%%%% LOAD THE SAME TRAJECTORIES AS THE TRIPTYCH
n_conditions = numel(condition_names);
traj_data = cell(n_conditions, 1);
traj_time = cell(n_conditions, 1);
track_lengths = cell(n_conditions, 1);
used_files = cell(n_conditions, 1);
color_indices = cell(n_conditions, 1);
n_color_slots = zeros(n_conditions, 1);

for c = 1:n_conditions
    data_dir = fullfile(base_dir, condition_dirs(c));
    [traj_data{c}, traj_time{c}, track_lengths{c}, used_files{c}] = ...
        load_snapshot_trajectories(data_dir, pix_mm_per_px, dt, ...
        use_smoothing, go, gw, min_timesteps_keep);

    % Keep the original triptych color assignment after removing a file.
    n_color_slots(c) = numel(traj_data{c});
    color_indices{c} = (1:n_color_slots(c))';
    manual_exclusion = ismember(used_files{c}, excluded_files{c});
    repeat_exclusion = remove_repeat_trajectories & ...
        ismember(used_files{c}, repeat_files{c});
    if any(manual_exclusion)
        fprintf('Excluding rejected %s trajectory: %s\n', condition_names(c), ...
            strjoin(used_files{c}(manual_exclusion), ", "));
    end
    if any(repeat_exclusion)
        fprintf('Excluding repeated %s trajectory: %s\n', condition_names(c), ...
            strjoin(used_files{c}(repeat_exclusion), ", "));
    end
    keep = ~(manual_exclusion | repeat_exclusion);
    traj_data{c} = traj_data{c}(keep);
    traj_time{c} = traj_time{c}(keep);
    track_lengths{c} = track_lengths{c}(keep);
    used_files{c} = used_files{c}(keep);
    color_indices{c} = color_indices{c}(keep);
end

clear warning_cleanup;

for c = selected_conditions
    fprintf('%s trajectories kept: %d\n', condition_names(c), numel(traj_data{c}));
    fprintf('%s files: %s\n', condition_names(c), strjoin(used_files{c}, ", "));
end

%%%% FIXED SQUARE CROP
% Coordinates are displayed in millimeters. All panels use the same 0--1
% mm window, with grid lines on all four box edges.
all_x = [];
all_y = [];
for c = 1:n_conditions
    for k = 1:numel(traj_data{c})
        all_x = [all_x; traj_data{c}{k}(:,1)]; %#ok<AGROW>
        all_y = [all_y; traj_data{c}{k}(:,2)]; %#ok<AGROW>
    end
end

if isempty(all_x)
    error('No trajectories were available after filtering.');
end

axis_limits = crop_limits_mm;

%%%% OPTIONALLY SPREAD FINAL-FRAME TRAJECTORIES BY TRANSLATION
trajectory_translations = cell(n_conditions, 1);
if minimize_trajectory_overlap
    last_frame_index = round(max(snapshot_times_s) / dt) + 1;
    for c = selected_conditions
        [traj_data{c}, trajectory_translations{c}, overlap_before, overlap_after] = ...
            spread_trajectories_by_translation(traj_data{c}, last_frame_index, ...
            crop_limits_mm, overlap_clearance_mm, overlap_grid_step_mm, ...
            overlap_candidate_count, overlap_optimization_passes);
        fprintf(['%s overlap optimization: %d -> %d overlapping trajectory ', ...
            'pairs; %d -> %d shared grid cells.\n'], condition_names(c), ...
            overlap_before.pair_count, overlap_after.pair_count, ...
            overlap_before.shared_cells, overlap_after.shared_cells);
    end
else
    for c = selected_conditions
        trajectory_translations{c} = zeros(numel(traj_data{c}), 2);
    end
end

%%%% COLORS MATCH animation_pico_triptych.m
colors = cell(n_conditions, 1);
for c = 1:n_conditions
    colors{c} = lines(max(n_color_slots(c), 1));
end

%%%% LAYOUT
switch lower(condition_layout)
    case "rows"
        n_rows = numel(selected_conditions);
        n_cols = numel(snapshot_times_s);
        figure_size = [100, 100, 1650, 1100];
        tile_index = @(condition_index, time_index) ...
            (condition_index - 1) * n_cols + time_index;
    case "columns"
        n_rows = numel(snapshot_times_s);
        n_cols = numel(selected_conditions);
        figure_size = [100, 100, 1150, 1650];
        tile_index = @(condition_index, time_index) ...
            (time_index - 1) * n_cols + condition_index;
    otherwise
        error('Unsupported condition_layout: %s. Use "rows" or "columns".', condition_layout);
end

fig = figure('Color', 'w', 'Position', figure_size);
layout = tiledlayout(fig, n_rows, n_cols, ...
    'TileSpacing', 'none', 'Padding', 'compact');
layout.OuterPosition = [0.08, 0.065, 0.91, 0.87];
axes_grid = gobjects(numel(selected_conditions), numel(snapshot_times_s));

axis_font_size = 14;
axis_label_font_size = 17;
grid_label_font_size = 19;

for condition_index = 1:numel(selected_conditions)
    c = selected_conditions(condition_index);

    for time_index = 1:numel(snapshot_times_s)
        snapshot_time = snapshot_times_s(time_index);
        frame_index = round(snapshot_time / dt) + 1; % frame 1 is t = 0

        ax = nexttile(layout, tile_index(condition_index, time_index));
        axes_grid(condition_index, time_index) = ax;
        hold(ax, 'on');

        for k = 1:numel(traj_data{c})
            xy = traj_data{c}{k};
            n_points = size(xy, 1);
            last_visible = min(frame_index, n_points);
            trajectory_color = colors{c}(color_indices{c}(k),:);

            if show_trails
                if isinf(trail_length)
                    first_visible = 1;
                else
                    first_visible = max(1, last_visible - trail_length + 1);
                end

                plot(ax, xy(first_visible:last_visible,1), ...
                    xy(first_visible:last_visible,2), '-', ...
                    'LineWidth', 1.5, 'Color', trajectory_color);
            end

            % In the GIF, the head disappears after a trajectory ends while
            % its completed trail remains visible.
            if frame_index <= n_points
                plot(ax, xy(frame_index,1), xy(frame_index,2), 'o', ...
                    'MarkerSize', 6, ...
                    'MarkerFaceColor', trajectory_color, ...
                    'MarkerEdgeColor', trajectory_color);
            end
        end

        axis(ax, 'equal');
        pbaspect(ax, [1 1 1]);
        xlim(ax, axis_limits(1:2));
        ylim(ax, axis_limits(3:4));
        xticks(ax, axis_limits(1):major_tick_step_mm:axis_limits(2));
        yticks(ax, axis_limits(3):major_tick_step_mm:axis_limits(4));
        grid(ax, 'on');
        box(ax, 'on');
        ax.LineWidth = 1.5;
        ax.FontSize = axis_font_size;

        tile = tile_index(condition_index, time_index);
        tile_row = ceil(tile / n_cols);
        tile_col = mod(tile - 1, n_cols) + 1;
        if tile_row < n_rows
            ax.XTickLabel = [];
        end
        if tile_col > 1
            ax.YTickLabel = [];
        elseif n_rows > 1
            % Touching panels share a border but represent different y
            % endpoints. Suppress both labels at each seam so they do not
            % print on top of one another.
            y_tick_labels = arrayfun(@(value) sprintf('%g', value), ...
                ax.YTick, 'UniformOutput', false);
            if tile_row < n_rows
                y_tick_labels{1} = '';
            end
            if tile_row > 1
                y_tick_labels{end} = '';
            end
            ax.YTickLabel = y_tick_labels;
        end
    end
end

xlabel(layout, 'x (mm)', 'FontSize', axis_label_font_size, 'FontWeight', 'bold');
ylabel(layout, 'y (mm)', 'FontSize', axis_label_font_size, 'FontWeight', 'bold');

drawnow;
add_grid_labels(fig, axes_grid, condition_layout, condition_names(selected_conditions), ...
    snapshot_times_s, grid_label_font_size);

exportgraphics(fig, output_png, ...
    'Resolution', output_resolution_dpi, 'BackgroundColor', 'white');
fprintf('Snapshot grid saved as: %s\n', output_png);

function [traj_data, traj_time, track_lengths, used_files] = ...
    load_snapshot_trajectories(data_dir, pix_mm_per_px, dt, ...
    use_smoothing, go, gw, min_timesteps_keep)

    if ~isfolder(data_dir)
        error('Input folder not found: %s', data_dir);
    end

    files = dir(fullfile(data_dir, '*.csv'));
    files = files(~[files.isdir]);
    if isempty(files)
        error('No CSV files found in %s.', data_dir);
    end

    file_numbers = nan(numel(files), 1);
    for k = 1:numel(files)
        [~, name, ~] = fileparts(files(k).name);
        file_numbers(k) = str2double(name);
    end
    [~, order] = sort(file_numbers);
    files = files(order);

    traj_data = {};
    traj_time = {};
    track_lengths = [];
    used_files = strings(0, 1);

    for k = 1:numel(files)
        filename = fullfile(files(k).folder, files(k).name);

        if files(k).bytes == 0
            fprintf('Skipping empty file: %s\n', filename);
            continue;
        end

        try
            table_data = readtable(filename);
        catch read_error
            fprintf('Skipping unreadable file: %s (%s)\n', filename, read_error.message);
            continue;
        end

        variables = string(table_data.Properties.VariableNames);
        if ~all(ismember(["X", "Y"], variables))
            fprintf('Skipping file with missing X/Y columns: %s\n', filename);
            continue;
        end

        x_px = table_data.X(:);
        y_px = table_data.Y(:);
        if ismember("FrameNumber", variables)
            frame_number = table_data.FrameNumber(:);
        else
            frame_number = (0:numel(x_px)-1)';
        end

        valid = isfinite(frame_number) & isfinite(x_px) & isfinite(y_px);
        frame_number = frame_number(valid);
        x_px = x_px(valid);
        y_px = y_px(valid);

        if numel(x_px) < min_timesteps_keep
            fprintf('Discarding short trajectory: %s (%d timesteps)\n', ...
                filename, numel(x_px));
            continue;
        end

        [frame_number, unique_index] = unique(frame_number, 'stable');
        x_px = x_px(unique_index);
        y_px = y_px(unique_index);

        x = x_px * pix_mm_per_px;
        y = y_px * pix_mm_per_px;
        t = (frame_number - frame_number(1)) * dt;

        if use_smoothing && numel(x) >= gw
            x = sgolayfilt(x, go, gw);
            y = sgolayfilt(y, go, gw);
        end

        if any(~isfinite(x)) || any(~isfinite(y)) || any(~isfinite(t))
            fprintf('Skipping file with non-finite values: %s\n', filename);
            continue;
        end

        traj_data{end+1,1} = [x(:), y(:)]; %#ok<AGROW>
        traj_time{end+1,1} = t(:); %#ok<AGROW>
        track_lengths(end+1,1) = numel(x); %#ok<AGROW>
        used_files(end+1,1) = string(files(k).name); %#ok<AGROW>
    end
end

function [translated_data, translations_mm, stats_before, stats_after] = ...
    spread_trajectories_by_translation(traj_data, last_frame_index, ...
    crop_limits, clearance_mm, grid_step_mm, candidate_count, n_passes)

    n_trajectories = numel(traj_data);
    translated_data = traj_data;
    translations_mm = zeros(n_trajectories, 2);
    if n_trajectories == 0
        stats_before = struct('pair_count', 0, 'shared_cells', 0);
        stats_after = stats_before;
        return;
    end

    if crop_limits(2) <= crop_limits(1) || crop_limits(4) <= crop_limits(3)
        error('crop_limits_mm must have increasing x and y bounds.');
    end
    if clearance_mm < 0 || grid_step_mm <= 0
        error('Overlap clearance must be nonnegative and grid step must be positive.');
    end

    n_grid_x = round((crop_limits(2) - crop_limits(1)) / grid_step_mm) + 1;
    n_grid_y = round((crop_limits(4) - crop_limits(3)) / grid_step_mm) + 1;
    base_pixels = cell(n_trajectories, 1);
    shift_bounds = zeros(n_trajectories, 4); % [sx_min sx_max sy_min sy_max]
    path_lengths = zeros(n_trajectories, 1);
    buffer_radius = 0.5 * clearance_mm;

    for k = 1:n_trajectories
        xy = traj_data{k};
        final_path = xy(1:min(last_frame_index, size(xy, 1)), :);
        path_lengths(k) = sum(vecnorm(diff(final_path, 1, 1), 2, 2));

        % Keep the buffered centerline inside the crop when possible. If a
        % trajectory is too large for that margin, require only its actual
        % centerline to remain inside.
        x_bounds = integer_shift_bounds(final_path(:,1), crop_limits(1:2), ...
            buffer_radius, grid_step_mm);
        y_bounds = integer_shift_bounds(final_path(:,2), crop_limits(3:4), ...
            buffer_radius, grid_step_mm);
        if isempty(x_bounds)
            x_bounds = integer_shift_bounds(final_path(:,1), crop_limits(1:2), ...
                0, grid_step_mm);
        end
        if isempty(y_bounds)
            y_bounds = integer_shift_bounds(final_path(:,2), crop_limits(3:4), ...
                0, grid_step_mm);
        end
        if isempty(x_bounds) || isempty(y_bounds)
            error(['Trajectory %d cannot be translated wholly inside the ', ...
                'requested crop at the final snapshot.'], k);
        end
        shift_bounds(k,:) = [x_bounds, y_bounds];

        sampled_path = resample_polyline(final_path, 0.5 * grid_step_mm);
        base_pixels{k} = buffered_path_pixels(sampled_path, crop_limits, ...
            grid_step_mm, buffer_radius);
    end

    zero_shifts = zeros(n_trajectories, 2);
    original_masks = build_shifted_masks(base_pixels, zero_shifts, ...
        n_grid_y, n_grid_x);
    stats_before = trajectory_overlap_stats(original_masks, n_grid_y * n_grid_x);

    % Place long trajectories first, then refine every placement while all
    % other paths are held fixed. Ties favor the smallest displacement from
    % the measured initial condition.
    [~, placement_order] = sort(path_lengths, 'descend');
    occupancy = zeros(n_grid_y * n_grid_x, 1, 'uint16');
    shifts = zeros(n_trajectories, 2);
    current_masks = cell(n_trajectories, 1);

    for order_index = 1:n_trajectories
        k = placement_order(order_index);
        preferred_shift = [clamp_integer(0, shift_bounds(k,1), shift_bounds(k,2)), ...
            clamp_integer(0, shift_bounds(k,3), shift_bounds(k,4))];
        [best_shift, best_mask] = best_translation_for_path(base_pixels{k}, ...
            shift_bounds(k,:), preferred_shift, occupancy, grid_step_mm, ...
            candidate_count, n_grid_y, n_grid_x);
        shifts(k,:) = best_shift;
        current_masks{k} = best_mask;
        occupancy(best_mask) = occupancy(best_mask) + 1;
    end

    for pass = 1:n_passes
        changed = false;
        if mod(pass, 2) == 0
            refinement_order = flipud(placement_order);
        else
            refinement_order = placement_order;
        end

        for order_index = 1:n_trajectories
            k = refinement_order(order_index);
            old_shift = shifts(k,:);
            old_mask = current_masks{k};
            occupancy(old_mask) = occupancy(old_mask) - 1;

            [best_shift, best_mask] = best_translation_for_path(base_pixels{k}, ...
                shift_bounds(k,:), old_shift, occupancy, grid_step_mm, ...
                candidate_count, n_grid_y, n_grid_x);
            shifts(k,:) = best_shift;
            current_masks{k} = best_mask;
            occupancy(best_mask) = occupancy(best_mask) + 1;
            changed = changed || any(best_shift ~= old_shift);
        end

        if ~changed
            break;
        end
    end

    translations_mm = shifts * grid_step_mm;
    for k = 1:n_trajectories
        translated_data{k} = traj_data{k} + translations_mm(k,:);
    end
    stats_after = trajectory_overlap_stats(current_masks, n_grid_y * n_grid_x);
end

function bounds = integer_shift_bounds(values, crop_bounds, margin, grid_step)
    lower = ceil((crop_bounds(1) + margin - min(values)) / grid_step - 1e-10);
    upper = floor((crop_bounds(2) - margin - max(values)) / grid_step + 1e-10);
    if lower > upper
        bounds = [];
    else
        bounds = [lower, upper];
    end
end

function value = clamp_integer(value, lower, upper)
    value = min(max(round(value), lower), upper);
end

function sampled = resample_polyline(xy, maximum_step)
    if size(xy, 1) < 2 || maximum_step <= 0
        sampled = xy;
        return;
    end

    sampled = zeros(0, 2);
    for point_index = 1:size(xy, 1)-1
        segment = xy(point_index+1,:) - xy(point_index,:);
        n_subsegments = max(1, ceil(norm(segment) / maximum_step));
        fractions = (0:n_subsegments-1)' / n_subsegments;
        sampled = [sampled; xy(point_index,:) + fractions .* segment]; %#ok<AGROW>
    end
    sampled = [sampled; xy(end,:)];
end

function pixels = buffered_path_pixels(xy, crop_limits, grid_step, radius_mm)
    columns = round((xy(:,1) - crop_limits(1)) / grid_step) + 1;
    rows = round((xy(:,2) - crop_limits(3)) / grid_step) + 1;
    radius_cells = ceil(radius_mm / grid_step);
    [column_offsets, row_offsets] = meshgrid(-radius_cells:radius_cells);
    disk = hypot(column_offsets, row_offsets) <= radius_mm / grid_step + 1e-10;
    column_offsets = column_offsets(disk)';
    row_offsets = row_offsets(disk)';

    expanded_rows = rows + row_offsets;
    expanded_columns = columns + column_offsets;
    pixels = unique([expanded_rows(:), expanded_columns(:)], 'rows');
end

function masks = build_shifted_masks(base_pixels, shifts, n_rows, n_columns)
    masks = cell(numel(base_pixels), 1);
    for k = 1:numel(base_pixels)
        masks{k} = shift_path_mask(base_pixels{k}, shifts(k,:), n_rows, n_columns);
    end
end

function mask = shift_path_mask(base_pixels, shift, n_rows, n_columns)
    rows = base_pixels(:,1) + shift(2);
    columns = base_pixels(:,2) + shift(1);
    inside = rows >= 1 & rows <= n_rows & columns >= 1 & columns <= n_columns;
    mask = unique(sub2ind([n_rows, n_columns], rows(inside), columns(inside)));
end

function values = candidate_integer_shifts(lower, upper, count, preferred)
    n_values = upper - lower + 1;
    if n_values <= count
        values = lower:upper;
    else
        values = unique(round(linspace(lower, upper, count)));
    end
    local_values = preferred + (-2:2);
    local_values = local_values(local_values >= lower & local_values <= upper);
    values = unique([values, local_values, clamp_integer(0, lower, upper)]);
end

function [best_shift, best_mask] = best_translation_for_path(base_pixels, ...
    bounds, preferred_shift, occupancy, grid_step, candidate_count, ...
    n_rows, n_columns)

    x_candidates = candidate_integer_shifts(bounds(1), bounds(2), ...
        candidate_count, preferred_shift(1));
    y_candidates = candidate_integer_shifts(bounds(3), bounds(4), ...
        candidate_count, preferred_shift(2));
    best_overlap = inf;
    best_motion = inf;
    best_shift = preferred_shift;
    best_mask = shift_path_mask(base_pixels, best_shift, n_rows, n_columns);

    for y_shift = y_candidates
        for x_shift = x_candidates
            candidate_shift = [x_shift, y_shift];
            candidate_mask = shift_path_mask(base_pixels, candidate_shift, ...
                n_rows, n_columns);
            overlap = sum(double(occupancy(candidate_mask)));
            motion = sum((candidate_shift * grid_step).^2);
            if overlap < best_overlap || ...
                    (overlap == best_overlap && motion < best_motion - 1e-14)
                best_overlap = overlap;
                best_motion = motion;
                best_shift = candidate_shift;
                best_mask = candidate_mask;
            end
        end
    end
end

function stats = trajectory_overlap_stats(masks, n_grid_cells)
    occupancy = zeros(n_grid_cells, 1, 'uint16');
    for k = 1:numel(masks)
        occupancy(masks{k}) = occupancy(masks{k}) + 1;
    end
    stats.shared_cells = nnz(occupancy > 1);
    stats.pair_count = 0;
    for first = 1:numel(masks)-1
        for second = first+1:numel(masks)
            if any(ismember(masks{first}, masks{second}))
                stats.pair_count = stats.pair_count + 1;
            end
        end
    end
end

function add_grid_labels(fig, axes_grid, condition_layout, ...
    selected_condition_names, snapshot_times_s, font_size)

    %#ok<INUSD> % fig is retained in the interface for future annotation-based labels
    label_args = {'Units', 'normalized', 'Clipping', 'off', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontSize', font_size, 'FontWeight', 'bold'};

    switch lower(condition_layout)
        case "rows"
            for time_index = 1:numel(snapshot_times_s)
                ax = axes_grid(1, time_index);
                text(ax, 0.5, 1.075, ...
                    sprintf('t = %g s', snapshot_times_s(time_index)), ...
                    label_args{:});
            end

            for condition_index = 1:numel(selected_condition_names)
                ax = axes_grid(condition_index, 1);
                text(ax, -0.145, 0.5, ...
                    char(selected_condition_names(condition_index)), ...
                    label_args{:}, 'Rotation', 90);
            end

        case "columns"
            for condition_index = 1:numel(selected_condition_names)
                ax = axes_grid(condition_index, 1);
                text(ax, 0.5, 1.075, ...
                    char(selected_condition_names(condition_index)), ...
                    label_args{:});
            end

            for time_index = 1:numel(snapshot_times_s)
                ax = axes_grid(1, time_index);
                text(ax, -0.145, 0.5, ...
                    sprintf('t = %g s', snapshot_times_s(time_index)), ...
                    label_args{:}, 'Rotation', 90);
            end
    end
end
