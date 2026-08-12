function [traj_data, track_lengths, used_ids, native_dts, common_dt] = ...
    load_standardized_resampled_trajectories(dir_data, file_pattern, traj_ids, ...
    use_smoothing, go, gw, common_dt, resample_method)
% Load standardized trajectories, smooth each at its native sampling, then
% resample all trajectories onto a shared uniform timestep for ensemble MSD.

    if nargin < 8 || isempty(resample_method)
        resample_method = 'pchip';
    end

    raw_traj_data = {};
    native_dts = [];
    used_ids = [];

    warning('off', 'all');
    for n = traj_ids
        filename = fullfile(dir_data, sprintf(file_pattern, n));

        if ~isfile(filename)
            fprintf('Missing file, skipping: %s\n', filename);
            continue;
        end

        [x, y, t, dt_native, meta] = load_standardized_trajectory(filename);

        if ~meta.has_time || ~isfinite(dt_native) || dt_native <= 0
            fprintf('Missing valid Time_s/dt, skipping: %s\n', filename);
            continue;
        end

        if numel(x) < 2
            fprintf('Trajectory too short, skipping: %s\n', filename);
            continue;
        end

        if use_smoothing && numel(x) >= gw
            x = sgolayfilt(x, go, gw);
            y = sgolayfilt(y, go, gw);
        end

        raw_traj_data{end+1,1} = [t(:), x(:), y(:)];
        native_dts(end+1,1) = dt_native;
        used_ids(end+1,1) = n;
    end
    warning('on', 'all');

    if isempty(raw_traj_data)
        error('No trajectories were loaded.');
    end

    if nargin < 7 || isempty(common_dt)
        common_dt = max(native_dts);
    end

    if ~isfinite(common_dt) || common_dt <= 0
        error('Invalid common_dt: %g', common_dt);
    end

    traj_data = {};
    track_lengths = [];
    keep_mask = false(numel(raw_traj_data), 1);

    for k = 1:numel(raw_traj_data)
        txy = raw_traj_data{k};
        t = txy(:,1);
        x = txy(:,2);
        y = txy(:,3);

        duration = t(end) - t(1);
        n_uniform = floor(duration / common_dt) + 1;

        if n_uniform < 2
            fprintf('Trajectory %d too short after resampling to dt=%.6g s, skipping.\n', used_ids(k), common_dt);
            continue;
        end

        t_uniform = t(1) + (0:n_uniform-1)' * common_dt;
        x_uniform = interp1(t, x, t_uniform, resample_method);
        y_uniform = interp1(t, y, t_uniform, resample_method);

        if any(~isfinite(x_uniform)) || any(~isfinite(y_uniform))
            fprintf('Interpolation failed, skipping trajectory %d.\n', used_ids(k));
            continue;
        end

        traj_data{end+1,1} = [x_uniform(:), y_uniform(:)];
        track_lengths(end+1,1) = numel(t_uniform);
        keep_mask(k) = true;
    end

    used_ids = used_ids(keep_mask);
    native_dts = native_dts(keep_mask);

    if isempty(traj_data)
        error('No trajectories remained after resampling.');
    end
end
