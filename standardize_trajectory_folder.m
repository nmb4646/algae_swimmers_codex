function standardize_trajectory_folder(input_dir, output_dir, pix_m_per_px, native_dt_s, target_dt_s, varargin)
% Standardize every CSV trajectory in a folder and optionally resample to a
% shared timestep before writing into the standardized-data folder.
%
% Required inputs:
%   input_dir     - folder containing raw CSV trajectories
%   output_dir    - folder to write standardized CSVs into
%   pix_m_per_px  - pixel size in meters per pixel
%   native_dt_s   - original timestep in seconds
%   target_dt_s   - output timestep in seconds
%
% Name-value options:
%   'FilePattern'   : default '*.csv'
%   'ZeroOrigin'    : default true
%   'InterpMethod'  : default 'pchip'
%   'SgolayOrder'   : default 3
%   'SgolayWindow'  : default 11
%   'DoSmoothing'   : default false

    p = inputParser;
    addRequired(p, 'input_dir', @(x) ischar(x) || isstring(x));
    addRequired(p, 'output_dir', @(x) ischar(x) || isstring(x));
    addRequired(p, 'pix_m_per_px', @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addRequired(p, 'native_dt_s', @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addRequired(p, 'target_dt_s', @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(p, 'FilePattern', '*.csv', @(x) ischar(x) || isstring(x));
    addParameter(p, 'ZeroOrigin', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'InterpMethod', 'pchip', @(x) ischar(x) || isstring(x));
    addParameter(p, 'SgolayOrder', 3, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'SgolayWindow', 11, @(x) isnumeric(x) && isscalar(x) && x >= 3);
    addParameter(p, 'DoSmoothing', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
    parse(p, input_dir, output_dir, pix_m_per_px, native_dt_s, target_dt_s, varargin{:});

    file_pattern = string(p.Results.FilePattern);
    zero_origin = logical(p.Results.ZeroOrigin);
    interp_method = char(string(p.Results.InterpMethod));
    sgolay_order = p.Results.SgolayOrder;
    sgolay_window = p.Results.SgolayWindow;
    do_smoothing = logical(p.Results.DoSmoothing);

    input_dir = char(string(input_dir));
    output_dir = char(string(output_dir));

    if ~isfolder(input_dir)
        error('Input folder not found: %s', input_dir);
    end
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end

    files = dir(fullfile(input_dir, char(file_pattern)));
    files = files(~[files.isdir]);
    if isempty(files)
        error('No files matching %s found in %s', file_pattern, input_dir);
    end

    summary_rows = cell(numel(files), 7);

    for k = 1:numel(files)
        input_path = fullfile(files(k).folder, files(k).name);
        T = readtable(input_path);

        vars = string(T.Properties.VariableNames);
        if ~all(ismember(["X", "Y"], vars))
            error('File %s is missing X/Y columns.', input_path);
        end

        x_px = T.X(:);
        y_px = T.Y(:);
        n_points_raw = numel(x_px);

        if ismember("FrameNumber", vars)
            frame_number_native = T.FrameNumber(:);
        else
            frame_number_native = (0:n_points_raw-1)';
        end

        valid_rows = isfinite(frame_number_native) & isfinite(x_px) & isfinite(y_px);
        n_dropped = sum(~valid_rows);
        if n_dropped > 0
            fprintf('Dropping %d non-finite rows from %s\n', n_dropped, files(k).name);
        end

        frame_number_native = frame_number_native(valid_rows);
        x_px = x_px(valid_rows);
        y_px = y_px(valid_rows);

        if numel(x_px) < 2
            fprintf('Skipping %s: fewer than 2 valid rows after filtering\n', files(k).name);
            continue;
        end

        [frame_number_native, unique_idx] = unique(frame_number_native, 'stable');
        x_px = x_px(unique_idx);
        y_px = y_px(unique_idx);

        n_points = numel(x_px);

        if zero_origin
            x_px = x_px - x_px(1);
            y_px = y_px - y_px(1);
        end

        x_m = x_px * pix_m_per_px;
        y_m = y_px * pix_m_per_px;
        t_native = (frame_number_native - frame_number_native(1)) * native_dt_s;

        if do_smoothing && n_points >= sgolay_window
            x_m = sgolayfilt(x_m, sgolay_order, sgolay_window);
            y_m = sgolayfilt(y_m, sgolay_order, sgolay_window);
        end

        if any(~isfinite(t_native)) || any(~isfinite(x_m)) || any(~isfinite(y_m))
            fprintf('Skipping %s: non-finite values remain after preprocessing\n', files(k).name);
            continue;
        end

        duration = t_native(end) - t_native(1);
        n_resampled = floor(duration / target_dt_s) + 1;
        if n_resampled < 2
            fprintf('Skipping %s: too short for target dt %.6g s\n', files(k).name, target_dt_s);
            continue;
        end

        t_out = t_native(1) + (0:n_resampled-1)' * target_dt_s;
        x_out = interp1(t_native, x_m, t_out, interp_method);
        y_out = interp1(t_native, y_m, t_out, interp_method);

        if any(~isfinite(x_out)) || any(~isfinite(y_out))
            fprintf('Skipping %s: interpolation produced non-finite values\n', files(k).name);
            continue;
        end

        output_path = fullfile(output_dir, files(k).name);
        frame_number_out = (0:n_resampled-1)';
        Tout = table(frame_number_out, t_out, x_out, y_out, ...
            'VariableNames', {'FrameNumber', 'Time_s', 'X', 'Y'});
        writetable(Tout, output_path);

        summary_rows(k, :) = {string(files(k).name), n_points, n_resampled, pix_m_per_px, native_dt_s, target_dt_s, string(output_path)};
        fprintf('Standardized %s -> %s\n', input_path, output_path);
    end

    summary_table = cell2table(summary_rows, 'VariableNames', ...
        {'file_name', 'n_input_points', 'n_output_points', 'pix_m_per_px', 'native_dt_s', 'target_dt_s', 'output_file'});
    writetable(summary_table, fullfile(output_dir, 'standardization_summary.csv'));

    fprintf('\nWrote %d standardized trajectories to %s\n', numel(files), output_dir);
    fprintf('Target dt: %.6g s\n', target_dt_s);
end
