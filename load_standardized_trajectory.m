function [x, y, t, dt, meta] = load_standardized_trajectory(filename)
% Load a standardized trajectory with coordinates in meters and time in seconds.

    T = readtable(filename);
    vars = string(T.Properties.VariableNames);

    if ~all(ismember(["X", "Y"], vars))
        error("File %s is missing X/Y columns.", filename);
    end

    x = T.X(:);
    y = T.Y(:);

    if ismember("Time_s", vars)
        t = T.Time_s(:);
        if numel(t) >= 2
            dt = median(diff(t), "omitnan");
        else
            dt = NaN;
        end
    else
        t = NaN(size(x));
        dt = NaN;
    end

    meta = struct();
    meta.filename = filename;
    meta.n_points = numel(x);
    meta.has_time = ismember("Time_s", vars);
    meta.variables = vars;
end
