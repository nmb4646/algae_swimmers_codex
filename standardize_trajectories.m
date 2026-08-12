clear; clc;

manifest_path = "./config/trajectory_standardization_manifest.csv";
output_root = "./data_standardized";
default_zero_origin = true;

if ~isfile(manifest_path)
    error("Manifest not found: %s", manifest_path);
end

if ~isfolder(output_root)
    mkdir(output_root);
end

opts = detectImportOptions(manifest_path, "TextType", "string");
M = readtable(manifest_path, opts);

required_cols = ["input_file", "pix_m_per_px", "dt_s"];
missing_cols = required_cols(~ismember(required_cols, string(M.Properties.VariableNames)));
if ~isempty(missing_cols)
    error("Manifest is missing required columns: %s", strjoin(missing_cols, ", "));
end

if ~ismember("output_subdir", string(M.Properties.VariableNames))
    M.output_subdir = strings(height(M), 1);
end
if ~ismember("output_file", string(M.Properties.VariableNames))
    M.output_file = strings(height(M), 1);
end
if ~ismember("zero_origin", string(M.Properties.VariableNames))
    M.zero_origin = repmat(default_zero_origin, height(M), 1);
end

summary_rows = cell(height(M), 10);

for i = 1:height(M)
    input_file = strtrim(M.input_file(i));
    output_subdir = strtrim(M.output_subdir(i));
    output_file = strtrim(M.output_file(i));
    pix_m_per_px = M.pix_m_per_px(i);
    dt_s = M.dt_s(i);
    zero_origin = parse_logical_flag(M.zero_origin(i), default_zero_origin);

    if input_file == ""
        error("Row %d has an empty input_file.", i);
    end
    if ~isfile(input_file)
        error("Input file not found on row %d: %s", i, input_file);
    end
    if ~isfinite(pix_m_per_px) || pix_m_per_px <= 0
        error("Row %d has invalid pix_m_per_px: %g", i, pix_m_per_px);
    end
    if ~isfinite(dt_s) || dt_s <= 0
        error("Row %d has invalid dt_s: %g", i, dt_s);
    end

    T = readtable(input_file);
    if ~all(ismember(["X", "Y"], string(T.Properties.VariableNames)))
        error("Input file is missing X/Y columns: %s", input_file);
    end

    x_px = T.X(:);
    y_px = T.Y(:);
    n_points = numel(x_px);

    if ismember("FrameNumber", string(T.Properties.VariableNames))
        frame_number = T.FrameNumber(:);
    else
        frame_number = (0:n_points-1)';
    end

    if zero_origin
        x_px = x_px - x_px(1);
        y_px = y_px - y_px(1);
    end

    x_m = x_px * pix_m_per_px;
    y_m = y_px * pix_m_per_px;
    time_s = (frame_number - frame_number(1)) * dt_s;

    if output_file == ""
        [~, base_name, ~] = fileparts(input_file);
        output_file = base_name + "_standardized.csv";
    end

    if output_subdir == ""
        out_dir = output_root;
    else
        out_dir = fullfile(output_root, output_subdir);
    end

    if ~isfolder(out_dir)
        mkdir(out_dir);
    end

    output_path = fullfile(out_dir, output_file);

    Tout = table(frame_number, time_s, x_m, y_m, ...
        'VariableNames', {'FrameNumber', 'Time_s', 'X', 'Y'});
    writetable(Tout, output_path);

    summary_rows(i, :) = { ...
        string(input_file), ...
        string(output_path), ...
        n_points, ...
        pix_m_per_px, ...
        dt_s, ...
        zero_origin, ...
        min(time_s), ...
        max(time_s), ...
        min(x_m.^2 + y_m.^2), ...
        max(x_m.^2 + y_m.^2)};

    fprintf("Standardized %s -> %s\n", input_file, output_path);
end

summary_table = cell2table(summary_rows, 'VariableNames', { ...
    'input_file', 'output_file', 'n_points', 'pix_m_per_px', 'dt_s', ...
    'zero_origin', 'time_start_s', 'time_end_s', 'min_radius_sq_m2', 'max_radius_sq_m2'});

writetable(summary_table, fullfile(output_root, "standardization_summary.csv"));
fprintf("\nWrote summary: %s\n", fullfile(output_root, "standardization_summary.csv"));

function tf = parse_logical_flag(value, default_value)
    if islogical(value)
        tf = value;
        return;
    end

    if isnumeric(value)
        tf = logical(value);
        return;
    end

    value_str = lower(strtrim(string(value)));
    if value_str == ""
        tf = default_value;
    elseif any(value_str == ["true", "t", "yes", "y", "1"])
        tf = true;
    elseif any(value_str == ["false", "f", "no", "n", "0"])
        tf = false;
    else
        error("Could not parse logical value: %s", value_str);
    end
end
