clear; clc; close all;

%%%% FIGURE 2j-m: PICOALGAE MOTION IN 0% AND 2% MUCIN
% This script builds all four panels from one curated set of independent
% swimmers. Coordinates are loaded once, then separately preprocessed for
% the MSD and geometry analyses. The central MSD calculation intentionally
% matches total_ensemble_pico_mucin_comparator.m (pooled geometric mean).
% j: geometric MSD; k: curvature PDF; l: directional autocorrelation;
% m: PDF of modeled long-time effective diffusion, one value per swimmer.
% Choose either a direct fit of the paper's algebraic single-swimmer MSD
% (Eq. 5) or the earlier speed/turn-rate/DACF route. Both extrapolate a
% short observed track; panel m is not a measured long-time MSD slope.
%
% Circle-swimmer model: dtheta = omega*dt + sqrt(2*Dr)*dW,
% C(tau) = exp(-Dr*tau)*cos(omega*tau), and in TWO dimensions
% D_eff = D_t + v^2*Dr / (2*(Dr^2 + omega^2)).
% Marine et al., Phys. Rev. E 87, 052305 (2013), Eq. (6),
% https://doi.org/10.1103/PhysRevE.87.052305
% The direct MSD route fits their Eq. (5), equivalently
% MSD(t) = 4*D_t*t + 2*v^2*Re[(exp(-z*t)-1+z*t)/z^2], z=D_r-1i*omega.
%
% Panel m uses a Gaussian KDE of z=log(D_eff/(1 um^2/s)), transformed back
% as p_D(D)=p_z(log(D))/D. One common log-bandwidth is selected by balanced
% within-condition leave-one-swimmer-out predictive likelihood. This is a
% density of D, NOT of log(D), and is independent of the axis-scale option.
% Jones et al., Positive data KDE, https://arxiv.org/abs/1804.08365

%%%% DATA AND CURATION SETTINGS
base_dir = './data/Picoalgae in 0-2% mucin/';
condition_dirs = ["picoalgae 0% mucin"; "picoalgae 2% mucin"];
condition_names = ["0% mucin"; "2% mucin"];
condition_colors = [0.12, 0.47, 0.71; ...   % blue
                    1.00, 0.25, 0.39];      % red/pink

fps = 15;
dt = 1 / fps;
pixel_size_um = 0.48;

% Require at least 4 s of data per swimmer; longer MSD lags use only the
% tracks that actually reach those lags (see msd_track_count in results).
analysis_duration_s = 4.0;

% These are the independently identified outliers and repeated tracks from
% the original 0--2% delivery. The retained member of each repeat group is
% not listed here.
manual_exclusions = { ...
    ["7.csv"; "11.csv"]; ...
    strings(0,1)};

repeat_exclusions = { ...
    ["5.csv"; "9.csv"; "10.csv"; "15.csv"; "16.csv"; ...
     "27.csv"; "37.csv"; "41.csv"; "42.csv"; "44.csv"]; ...
    ["14.csv"; "19.csv"; "23.csv"; "30.csv"]};

%%%% PREPROCESSING SETTINGS
% Preserve the comparator's MSD preprocessing exactly.
msd_sgolay_order = 3;
msd_smoothing_window_frames = 11;

% Geometry needs lighter smoothing: 5 frames = 0.333 s at 15 fps.
geometry_sgolay_order = 2;
geometry_smoothing_window_frames = 5;

%%%% PANEL-SPECIFIC SETTINGS
msd_max_lag_s = 25.0;
msd_loglog = true;                 % true = logarithmic time and MSD axes

curvature_max_um_inv = 0.50;
curvature_plot_max_um_inv = 0.30;
curvature_grid_points = 240;

dacf_max_lag_s = 4.0;
minimum_direction_step_um = 0.05;
show_dacf_fit = false;

% Use each swimmer's mean SIGNED turning rate as the constant omega;
% hold it fixed when fitting that swimmer's DACF for Dr. Averaging absolute
% instantaneous turning/curvature would count fluctuations as steady spin.
% Speed and heading use the same lightly smoothed steps as the DACF. All
% time origins in each complete track contribute; this is not a 4-s snippet.
diffusion_estimation_method = "single_particle_msd_fit"; % or "dacf_parameters"
msd_fit_use_sg3_11 = true; % match the standalone method-c PDF; panel j unchanged
msd_fit_max_lag_s = 25;
msd_fit_max_track_fraction = 0.5; % avoid very sparse long-lag estimates
msd_fit_min_lag_span_s = 4;
msd_fit_min_r_squared = 0;
msd_fit_warn_r_squared = 0.8;
diffusion_translational_baseline_um2_s = [0; 0]; % one D_t per condition
% D_t defaults to zero: no cell radius or mucin viscosity is being assumed.
% Supply independently justified values here to include passive diffusion.
diffusion_max_Dr_s_inv = 5 / dt;   % numerical search bound, NOT a noise prior
diffusion_min_fit_r_squared = 0;  % omit only fits worse than a constant mean
diffusion_warn_fit_r_squared = 0.50; % flag weak fits, do not remove them
diffusion_grid_points = 1200;    % full-support logarithmic integration grid
diffusion_bandwidth_log = [];   % [] = shared, cross-validated log-bandwidth
diffusion_cv_bandwidth_limits_log = []; % [] = data-scaled search bounds
diffusion_cv_grid_points = 161;
diffusion_bandwidth_sensitivity_factors = [0.75, 1, 1.25];
diffusion_bootstrap_samples = 1000;
diffusion_bootstrap_seed = 20260928;
diffusion_plot_xlim_um2_s = [0, 40]; % peak region; KDE remains normalized
diffusion_log_x = false;             % true = logarithmic D_eff axis
show_diffusion_peaks = true;         % modes of the full estimated PDFs, not means
show_diffusion_rug = false;          % cleaner manuscript panel
run_diffusion_self_tests = true;

%%%% OUTPUT SETTINGS
if diffusion_estimation_method == "single_particle_msd_fit"
    if msd_fit_use_sg3_11
        output_stem = 'analysis_2jklm_msd_fit_sg3_11_sem';
    else
        output_stem = 'analysis_2jklm_msd_fit';
    end
elseif diffusion_estimation_method == "dacf_parameters"
    output_stem = 'analysis_2jklm';
else
    error('Unknown diffusion estimation method: %s', diffusion_estimation_method);
end
output_resolution_dpi = 600;
save_individual_panels = true;
save_results_mat = true;
save_swimmer_parameters_csv = true;
save_diffusion_kde_diagnostics = true;
save_msd_fit_diagnostics = true;

if diffusion_estimation_method == "single_particle_msd_fit" && msd_fit_use_sg3_11
    % Match the shared bandwidth used by the existing SG3/11 a/b/c PDF.
    reference = load('compare_pico_diffusion_methods_abc_sg1_3_results.mat', ...
        'settings');
    diffusion_bandwidth_log = reference.settings.kde_bandwidth_log;
    assert(isfinite(diffusion_bandwidth_log) && diffusion_bandwidth_log > 0);
end

if run_diffusion_self_tests
    test_diffusion_model();
    test_circle_msd_fit();
end
if numel(diffusion_translational_baseline_um2_s) ~= numel(condition_names) || ...
        any(~isfinite(diffusion_translational_baseline_um2_s)) || ...
        any(diffusion_translational_baseline_um2_s < 0)
    error('Provide one finite, nonnegative D_t for each condition.');
end

%%%% LOAD ONE SHARED, CURATED POPULATION
n_conditions = numel(condition_names);
raw_trajectories = cell(n_conditions, 1);
used_files = cell(n_conditions, 1);
load_info = cell(n_conditions, 1);

for c = 1:n_conditions
    data_dir = fullfile(base_dir, condition_dirs(c));
    excluded = unique([manual_exclusions{c}; repeat_exclusions{c}]);

    [raw_trajectories{c}, used_files{c}, load_info{c}] = ...
        load_curated_trajectories(data_dir, pixel_size_um, dt, ...
        analysis_duration_s, excluded);

    fprintf('\n%s: retained %d independent swimmers\n', ...
        condition_names(c), numel(raw_trajectories{c}));
    fprintf('Files: %s\n', strjoin(used_files{c}, ', '));
end

%%%% PREPROCESS THE SAME SWIMMERS FOR EACH ANALYSIS FAMILY
msd_trajectories = cell(n_conditions, 1);
geometry_trajectories = cell(n_conditions, 1);

for c = 1:n_conditions
    msd_trajectories{c} = smooth_trajectory_set(raw_trajectories{c}, ...
        msd_sgolay_order, msd_smoothing_window_frames);
    geometry_trajectories{c} = smooth_trajectory_set(raw_trajectories{c}, ...
        geometry_sgolay_order, geometry_smoothing_window_frames);
end

%%%% COMPUTE PANEL J AND THE DACF FOR PANEL L / PER-SWIMMER MODELING
results = repmat(struct(), n_conditions, 1);

for c = 1:n_conditions
    results(c).condition = condition_names(c);
    results(c).files = used_files{c};
    results(c).load_info = load_info{c};

    [results(c).msd_time_s, results(c).msd_um2, ...
        results(c).msd_lower_um2, results(c).msd_upper_um2, ...
        results(c).msd_track_count] = ensemble_geometric_msd( ...
        msd_trajectories{c}, dt, msd_max_lag_s);

    [results(c).dacf_lag_s, results(c).dacf_mean, ...
        results(c).dacf_sem, results(c).dacf_matrix, ...
        results(c).dacf_pair_count_matrix] = ...
        ensemble_directional_autocorrelation(geometry_trajectories{c}, ...
        dt, dacf_max_lag_s, minimum_direction_step_um);

    results(c).dacf_fit = fit_damped_dacf(results(c).dacf_lag_s, ...
        results(c).dacf_mean, results(c).dacf_sem);

    fprintf(['%s DACF fit: D_r = %.3f s^-1, decay time = %.3f s, ', ...
        'Omega = %.3f rad/s, period = %.3f s, R^2 = %.3f\n'], ...
        condition_names(c), results(c).dacf_fit.Dr_s_inv, ...
        results(c).dacf_fit.decay_time_s, ...
        results(c).dacf_fit.omega_rad_s, ...
        results(c).dacf_fit.period_s, results(c).dacf_fit.r_squared);
end

%%%% COMPUTE PANEL K WITH ONE COMMON KDE BANDWIDTH
curvature_values = cell(n_conditions, 1);
all_curvature_values = [];

for c = 1:n_conditions
    curvature_values{c} = cell(numel(geometry_trajectories{c}), 1);
    for k = 1:numel(geometry_trajectories{c})
        kappa = menger_curvature(geometry_trajectories{c}{k});
        valid = isfinite(kappa) & kappa >= 0 & ...
            kappa <= curvature_max_um_inv;
        curvature_values{c}{k} = kappa(valid);
        all_curvature_values = [all_curvature_values; kappa(valid)]; %#ok<AGROW>
    end
end

if numel(all_curvature_values) < 2
    error('Too few valid curvature estimates were available for a PDF.');
end

curvature_grid = linspace(0, curvature_max_um_inv, ...
    curvature_grid_points);
% Estimate one bandwidth from all retained observations, then apply an
% explicit reflected KDE below. Doing the reflection manually preserves
% exactly straight (kappa = 0) observations, which MATLAB's built-in
% positive-support option rejects at the boundary.
[~, ~, curvature_bandwidth] = ksdensity(all_curvature_values, ...
    curvature_grid);

fprintf('Common curvature KDE bandwidth: %.5f um^-1\n', ...
    curvature_bandwidth);

for c = 1:n_conditions
    [results(c).curvature_pdf_mean, results(c).curvature_pdf_sem, ...
        results(c).curvature_pdf_matrix, ...
        results(c).curvature_valid_track_count] = ensemble_curvature_pdf( ...
        curvature_values{c}, curvature_grid, curvature_bandwidth);
    results(c).curvature_grid_um_inv = curvature_grid;
    results(c).curvature_peak_um_inv = peak_location_parabolic( ...
        curvature_grid, results(c).curvature_pdf_mean);

    fprintf('%s curvature peak: %.4f um^-1 (radius %.2f um)\n', ...
        condition_names(c), results(c).curvature_peak_um_inv, ...
        1 / results(c).curvature_peak_um_inv);
end

%%%% COMPUTE PANEL M: ONE LONG-TIME DIFFUSION COEFFICIENT PER SWIMMER
all_diffusion_values = [];
for c = 1:n_conditions
    if diffusion_estimation_method == "single_particle_msd_fit"
        fit_positions = raw_trajectories{c};
        if msd_fit_use_sg3_11
            fit_positions = msd_trajectories{c};
        end
        [parameters, fit_curves] = swimmer_msd_fit_parameters( ...
            fit_positions, used_files{c}, dt, ...
            msd_fit_max_lag_s, msd_fit_max_track_fraction, ...
            msd_fit_min_lag_span_s, ...
            diffusion_translational_baseline_um2_s(c), ...
            diffusion_max_Dr_s_inv, msd_fit_min_r_squared, ...
            msd_fit_warn_r_squared);
        results(c).diffusion_msd_fits = fit_curves;
    else
        [parameters, fitted_curves] = swimmer_diffusion_parameters( ...
            geometry_trajectories{c}, used_files{c}, curvature_values{c}, ...
            dt, minimum_direction_step_um, results(c).dacf_lag_s, ...
            results(c).dacf_matrix, results(c).dacf_pair_count_matrix, ...
            diffusion_translational_baseline_um2_s(c), ...
            diffusion_max_Dr_s_inv, diffusion_min_fit_r_squared, ...
            diffusion_warn_fit_r_squared);
        results(c).diffusion_dacf_fit_matrix = fitted_curves;
    end
    parameters = addvars(parameters, ...
        repmat(condition_names(c), height(parameters), 1), ...
        'Before', 1, 'NewVariableNames', 'condition');
    results(c).swimmer_parameters = parameters;
    included = parameters.included_in_diffusion_pdf;
    values = parameters.Deff_um2_s(included);
    if numel(values) < 2
        error('Too few model-compatible %s swimmers for a diffusion PDF.', ...
            condition_names(c));
    end
    results(c).diffusion_values_um2_s = values;
    results(c).diffusion_swimmer_count = numel(values);
    results(c).diffusion_mean_um2_s = mean(values);
    results(c).diffusion_median_um2_s = median(values);
    all_diffusion_values = [all_diffusion_values; values]; %#ok<AGROW>
    fprintf(['\n%s modeled diffusion: %d/%d swimmers; median %.2f, ', ...
        'IQR %.2f--%.2f um^2/s.\n'], condition_names(c), ...
        numel(values), height(parameters), median(values), ...
        prctile(values, 25), prctile(values, 75));
    fprintf('Excluded from diffusion PDF ONLY: %s\n', ...
        printable_list(parameters.file(~included)));
    if diffusion_estimation_method == "single_particle_msd_fit"
        fprintf('Included but weak MSD fit (R^2 < %.2f): %s\n', ...
            msd_fit_warn_r_squared, printable_list(parameters.file( ...
            included & parameters.fit_r_squared < msd_fit_warn_r_squared)));
        fprintf('MSD-fit decay time longer than fit window: %s\n', ...
            printable_list(parameters.file(included & ...
            ~parameters.noise_decay_observed)));
    else
        fprintf('Included but weak DACF fit (R^2 < %.2f): %s\n', ...
            diffusion_warn_fit_r_squared, printable_list(parameters.file( ...
            included & parameters.fit_r_squared < diffusion_warn_fit_r_squared)));
        fprintf('Decay time longer than record (extrapolation flag): %s\n', ...
            printable_list(parameters.file(included & ...
            ~parameters.noise_decay_observed)));
    end
end

%%%% POSITIVE-SUPPORT KDE WITH A SHARED, DATA-SELECTED BANDWIDTH
diffusion_groups = {results.diffusion_values_um2_s}';
[selected_bandwidth_log, diffusion_kde_cv] = select_log_kde_bandwidth( ...
    diffusion_groups, diffusion_cv_bandwidth_limits_log, diffusion_cv_grid_points);
if isempty(diffusion_bandwidth_log)
    diffusion_bandwidth_log = selected_bandwidth_log;
    diffusion_kde_cv.selection = "cross-validation";
else
    diffusion_kde_cv.selection = "manual override";
end
if ~isscalar(diffusion_bandwidth_log) || ~isfinite(diffusion_bandwidth_log) || ...
        diffusion_bandwidth_log <= 0
    error('The log-KDE bandwidth must be finite and positive.');
end
if any(~isfinite(diffusion_bandwidth_sensitivity_factors)) || ...
        any(diffusion_bandwidth_sensitivity_factors <= 0)
    error('Bandwidth sensitivity factors must be finite and positive.');
end
diffusion_kde_cv.used_bandwidth_log = diffusion_bandwidth_log;
% Cover six log-space standard deviations beyond all kernel centers, even
% for the broadest sensitivity curve. Do not normalize on a cropped range.
% A mixed log/linear grid resolves both the tails and the displayed peaks.
maximum_test_bandwidth = diffusion_bandwidth_log * ...
    max([1, diffusion_bandwidth_sensitivity_factors]);
minimum_test_bandwidth = diffusion_bandwidth_log * ...
    min([1, diffusion_bandwidth_sensitivity_factors]);
log_grid_limits = [min(log(all_diffusion_values))-6*maximum_test_bandwidth, ...
    max(log(all_diffusion_values))+6*maximum_test_bandwidth];
diffusion_axis_max = exp(log_grid_limits(2));
diffusion_grid_points = max(diffusion_grid_points, ...
    ceil(12*diff(log_grid_limits)/minimum_test_bandwidth)+1);
if isempty(diffusion_plot_xlim_um2_s)
    linear_grid_max = min(150, diffusion_axis_max);
else
    linear_grid_max = diffusion_plot_xlim_um2_s(end);
end
diffusion_grid = unique([0, exp(linspace(log_grid_limits(1), ...
    log_grid_limits(2), diffusion_grid_points)), ...
    linspace(0, linear_grid_max, 601)]);

for c = 1:n_conditions
    [results(c).diffusion_pdf_mean, results(c).diffusion_pdf_sem, ...
        results(c).diffusion_pdf_matrix] = swimmer_diffusion_pdf( ...
        results(c).diffusion_values_um2_s, diffusion_grid, ...
        diffusion_bandwidth_log);
    results(c).diffusion_grid_um2_s = diffusion_grid;
    [results(c).diffusion_pdf_lower, results(c).diffusion_pdf_upper, ...
        results(c).diffusion_bootstrap_peak_um2_s] = bootstrap_kde( ...
        results(c).diffusion_pdf_matrix, diffusion_grid, ...
        diffusion_bootstrap_samples, diffusion_bootstrap_seed+c-1);
    results(c).diffusion_peak_ci95_um2_s = prctile( ...
        results(c).diffusion_bootstrap_peak_um2_s, [2.5, 97.5]);
    results(c).diffusion_peak_um2_s = peak_location_parabolic( ...
        diffusion_grid, results(c).diffusion_pdf_mean);
    results(c).diffusion_peak_density = interp1(diffusion_grid, ...
        results(c).diffusion_pdf_mean, results(c).diffusion_peak_um2_s, 'pchip');
    results(c).diffusion_peak_at_boundary = ...
        results(c).diffusion_peak_um2_s == diffusion_grid(1) || ...
        results(c).diffusion_peak_um2_s == diffusion_grid(end);
    fprintf('%s log-KDE peak: %.2f um^2/s; bootstrap 95%% range %.2f--%.2f.\n', ...
        condition_names(c), results(c).diffusion_peak_um2_s, ...
        results(c).diffusion_peak_ci95_um2_s);
    n_factors = numel(diffusion_bandwidth_sensitivity_factors);
    results(c).diffusion_sensitivity_pdf = zeros(n_factors, numel(diffusion_grid));
    results(c).diffusion_sensitivity_peak_um2_s = nan(1, n_factors);
    for f = 1:n_factors
        density = swimmer_diffusion_pdf(results(c).diffusion_values_um2_s, ...
            diffusion_grid, diffusion_bandwidth_log * ...
            diffusion_bandwidth_sensitivity_factors(f));
        results(c).diffusion_sensitivity_pdf(f,:) = density;
        results(c).diffusion_sensitivity_peak_um2_s(f) = ...
            peak_location_parabolic(diffusion_grid, density);
    end
    results(c).diffusion_pdf_area = trapz(diffusion_grid, ...
        results(c).diffusion_pdf_mean);
    assert(abs(results(c).diffusion_pdf_area - 1) < 1e-3, ...
        'Diffusion PDF grid must resolve the normalized log-KDE.');
end
fprintf('Shared diffusion log-bandwidth: %.5f (%s).\n', ...
    diffusion_bandwidth_log, diffusion_kde_cv.selection);

%%%% MANUSCRIPT-STYLE FOUR-PANEL FIGURE
fig = figure('Color', 'w', 'Position', [50, 100, 1850, 470]);
layout = tiledlayout(fig, 1, 4, 'TileSpacing', 'compact', ...
    'Padding', 'loose');

% Preserve the paper's baseline limits, but never clip a calculated mean or
% its uncertainty band merely to imitate the draft artwork.
curvature_axis_max = 12;
for c = 1:n_conditions
    curvature_axis_max = max(curvature_axis_max, ...
        max(results(c).curvature_pdf_mean + results(c).curvature_pdf_sem, ...
        [], 'omitnan'));
end
curvature_axis_max = 2 * ceil(1.02 * curvature_axis_max / 2);

ax_j = nexttile(layout, 1); hold(ax_j, 'on');
plot_msd_panel(ax_j, results, condition_colors, condition_names, ...
    msd_loglog, dt, msd_max_lag_s);
style_panel(ax_j, 'j');

ax_k = nexttile(layout, 2); hold(ax_k, 'on');
for c = 1:n_conditions
    x = results(c).curvature_grid_um_inv;
    y = results(c).curvature_pdf_mean;
    sem = results(c).curvature_pdf_sem;
    visible = x <= curvature_plot_max_um_inv;
    x = x(visible);
    y = y(visible);
    sem = sem(visible);
    plot_band(ax_k, x, max(0, y - sem), y + sem, ...
        condition_colors(c,:), 0.18);
    plot(ax_k, x, y, 'Color', condition_colors(c,:), 'LineWidth', 2.0, ...
        'HandleVisibility', 'off');
    xline(ax_k, results(c).curvature_peak_um_inv, '--', ...
        'Color', condition_colors(c,:), 'LineWidth', 1.1, ...
        'HandleVisibility', 'off');
end
xlim(ax_k, [0, curvature_plot_max_um_inv]);
ylim(ax_k, [0, curvature_axis_max]);
xlabel(ax_k, 'Curvature \kappa (\mum^{-1})', 'Interpreter', 'tex');
ylabel(ax_k, 'Probability Density');
style_panel(ax_k, 'k');

ax_l = nexttile(layout, 3); hold(ax_l, 'on');
dacf_axis_values = 0;
for c = 1:n_conditions
    x = results(c).dacf_lag_s;
    y = results(c).dacf_mean;
    sem = results(c).dacf_sem;
    dacf_axis_values = [dacf_axis_values; y(:)-sem(:); ...
        y(:)+sem(:); y(:)]; %#ok<AGROW>
    plot_band(ax_l, x, y - sem, y + sem, condition_colors(c,:), 0.18);
    plot(ax_l, x, y, 'Color', condition_colors(c,:), 'LineWidth', 2.0, ...
        'HandleVisibility', 'off');
    plot_optional_dacf_fit(ax_l, x, ...
        results(c).dacf_fit.fitted_curve, condition_colors(c,:), ...
        show_dacf_fit);
    % Including the optional fit in the bounds also keeps the axes valid
    % if show_dacf_fit is switched on without changing any other setting.
    dacf_axis_values = [dacf_axis_values; ...
        results(c).dacf_fit.fitted_curve(:)]; %#ok<AGROW>
end
yline(ax_l, 0, ':', 'Color', [0.35, 0.35, 0.35], ...
    'LineWidth', 0.9, 'HandleVisibility', 'off');
xlim(ax_l, [0, dacf_max_lag_s]);
dacf_axis_values = dacf_axis_values(isfinite(dacf_axis_values));
dacf_padding = 0.05 * range(dacf_axis_values);
dacf_ymin = min(-0.5, floor(10 * (min(dacf_axis_values) - ...
    dacf_padding)) / 10);
dacf_ymax = max(1, ceil(10 * (max(dacf_axis_values) + ...
    dacf_padding)) / 10);
ylim(ax_l, [dacf_ymin, dacf_ymax]);
xlabel(ax_l, '\tau (s)', 'Interpreter', 'tex');
ylabel(ax_l, 'Direction Correlation');
style_panel(ax_l, 'l');

ax_m = nexttile(layout, 4); hold(ax_m, 'on');
if isempty(diffusion_plot_xlim_um2_s)
    diffusion_display_limits = [0, diffusion_axis_max];
else
    diffusion_display_limits = diffusion_plot_xlim_um2_s;
end
if numel(diffusion_display_limits) ~= 2 || ...
        any(~isfinite(diffusion_display_limits)) || ...
        diffusion_display_limits(1) < 0 || ...
        diffusion_display_limits(2) <= diffusion_display_limits(1)
    error('Diffusion plot limits must be [] or [nonnegative lower, larger upper].');
end
if diffusion_log_x
    if diffusion_display_limits(1) == 0 %#ok<UNRCH> enabled by diffusion_log_x
        % Zero cannot appear on a log axis. Pick a positive display bound
        % below every positive observed coefficient; keep the PDF in D
        % units (this is NOT a PDF of log(D)).
        positive_values = all_diffusion_values(all_diffusion_values > 0);
        diffusion_display_limits(1) = min([positive_values/2; ...
            diffusion_grid(find(diffusion_grid > 0, 1)); ...
            diffusion_display_limits(2)/10]);
    end
    set(ax_m, 'XScale', 'log');
end
diffusion_axis_values = [];
for c = 1:n_conditions
    x = results(c).diffusion_grid_um2_s;
    y = results(c).diffusion_pdf_mean;
    lower = max(0,y-results(c).diffusion_pdf_sem);
    upper = y+results(c).diffusion_pdf_sem;
    if diffusion_log_x
        positive = x > 0; %#ok<UNRCH> enabled by diffusion_log_x
        x = x(positive);
        y = y(positive);
        lower = lower(positive);
        upper = upper(positive);
    end
    plot_band(ax_m, x, max(0, lower), upper, ...
        condition_colors(c,:), 0.18);
    plot(ax_m, x, y, 'Color', condition_colors(c,:), 'LineWidth', 2.0, ...
        'HandleVisibility', 'off');
    diffusion_axis_values = [diffusion_axis_values; upper(:); y(:)]; %#ok<AGROW>
end
xlim(ax_m, diffusion_display_limits);
ylim(ax_m, padded_zero_ylim(diffusion_axis_values, 0.08));
xlabel(ax_m, 'Modeled D_{eff} (\mum^2/s)', 'Interpreter', 'tex');
ylabel(ax_m, 'Probability Density (s/\mum^2)', 'Interpreter', 'tex');
if show_diffusion_peaks
    for c = 1:n_conditions
        peak_D = results(c).diffusion_peak_um2_s;
        % Do not overwrite the y-axis with a vertical line at a boundary
        % mode, or try to draw a D=0 marker on a logarithmic axis.
        if peak_D > diffusion_display_limits(1) && ...
                peak_D < diffusion_display_limits(2)
            xline(ax_m, peak_D, '--', 'Color', condition_colors(c,:), ...
                'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
        peak_label = sprintf('Peak D_{eff} = %.1f', peak_D);
        if results(c).diffusion_peak_at_boundary
            peak_label = [peak_label, ' (boundary)']; %#ok<AGROW>
        end
        text(ax_m, 0.97, 0.86-0.08*(c-1), peak_label, ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'top', 'FontName', 'Arial', 'FontSize', 10, ...
            'Color', condition_colors(c,:));
    end
end
if show_diffusion_rug
    rug_scale = ax_m.YLim(2);
    for c = 1:n_conditions
        values = results(c).diffusion_values_um2_s;
        visible = values >= diffusion_display_limits(1) & ...
            values <= diffusion_display_limits(2);
        values = values(visible)';
        bottom = (c-1)*0.03*rug_scale;
        top = bottom+0.018*rug_scale;
        plot(ax_m, [values; values], ...
            [repmat(bottom, size(values)); repmat(top, size(values))], ...
            'Color', condition_colors(c,:), 'LineWidth', 1.0, ...
            'HandleVisibility', 'off');
    end
end
style_panel(ax_m, 'm');

%%%% EXPORT COMPOSITE, INDIVIDUAL PANELS, AND NUMERICAL RESULTS
drawnow;
savefig(fig, output_stem + ".fig");
if save_individual_panels
    exportgraphics(ax_j, output_stem + "_2j_msd.png", ...
        'Resolution', output_resolution_dpi);
    exportgraphics(ax_k, output_stem + "_2k_curvature.png", ...
        'Resolution', output_resolution_dpi);
    exportgraphics(ax_l, output_stem + "_2l_dacf.png", ...
        'Resolution', output_resolution_dpi);
    exportgraphics(ax_m, output_stem + "_2m_diffusion.png", ...
        'Resolution', output_resolution_dpi);
end
exportgraphics(fig, output_stem + ".png", ...
    'Resolution', output_resolution_dpi);
% Use a fresh graphics handle for the second full-figure export.
pdf_fig = openfig(output_stem + ".fig", 'invisible');
exportgraphics(pdf_fig, output_stem + ".pdf", 'ContentType', 'vector');
close(pdf_fig);

if save_swimmer_parameters_csv
    swimmer_parameters = vertcat(results.swimmer_parameters);
    writetable(swimmer_parameters, output_stem + "_swimmer_parameters.csv");
end

if save_diffusion_kde_diagnostics
    diagnostic_fig = plot_log_kde_diagnostics(results, diffusion_kde_cv, ...
        condition_colors, condition_names, diffusion_bandwidth_sensitivity_factors, ...
        diffusion_display_limits);
    exportgraphics(diagnostic_fig, output_stem + "_diffusion_kde_diagnostics.png", ...
        'Resolution', 300);
    exportgraphics(diagnostic_fig, output_stem + "_diffusion_kde_diagnostics.pdf", ...
        'ContentType', 'vector');

    intervals = vertcat(results.diffusion_peak_ci95_um2_s);
    kde_summary = table(condition_names, [results.diffusion_swimmer_count]', ...
        repmat(diffusion_bandwidth_log, n_conditions, 1), ...
        [results.diffusion_peak_um2_s]', intervals(:,1), intervals(:,2), ...
        [results.diffusion_mean_um2_s]', [results.diffusion_median_um2_s]', ...
        'VariableNames', {'condition', 'n_swimmers', 'bandwidth_log', ...
        'pdf_peak_um2_s', 'bootstrap_peak_p025_um2_s', ...
        'bootstrap_peak_p975_um2_s', 'mean_um2_s', 'median_um2_s'});
    sensitivity_peaks = vertcat(results.diffusion_sensitivity_peak_um2_s);
    for f = 1:numel(diffusion_bandwidth_sensitivity_factors)
        name = sprintf('peak_at_%03d_percent_bandwidth_um2_s', ...
            round(100*diffusion_bandwidth_sensitivity_factors(f)));
        kde_summary.(name) = sensitivity_peaks(:,f);
    end
    writetable(kde_summary, output_stem + "_diffusion_kde_summary.csv");
end

if save_msd_fit_diagnostics && ...
        diffusion_estimation_method == "single_particle_msd_fit"
    msd_fit_fig = plot_single_swimmer_msd_fit_diagnostics( ...
        results, condition_colors, condition_names);
    exportgraphics(msd_fit_fig, output_stem + "_fit_diagnostics.png", ...
        'Resolution', 300);
    exportgraphics(msd_fit_fig, output_stem + "_fit_diagnostics.pdf", ...
        'ContentType', 'vector');
end

if save_results_mat
    msd_fit_input = "raw curated positions; algebraic time-origin average";
    if msd_fit_use_sg3_11
        msd_fit_input = "SG3/11-smoothed positions; algebraic time-origin average";
    end
    if diffusion_estimation_method == "single_particle_msd_fit"
        omega_method = "joint algebraic single-swimmer MSD curve fit";
        fit_weights = "unweighted least squares across fitted MSD lags";
    else
        omega_method = "mean signed step-heading increment / dt";
        fit_weights = "valid direction-pair counts in swimmer DACF";
    end
    settings = struct( ...
        'base_dir', base_dir, ...
        'condition_dirs', condition_dirs, ...
        'condition_names', condition_names, ...
        'fps', fps, ...
        'dt_s', dt, ...
        'pixel_size_um', pixel_size_um, ...
        'analysis_duration_s', analysis_duration_s, ...
        'panel_order', ["MSD", "curvature PDF", "DACF", "modeled diffusion PDF"], ...
        'diffusion_estimation_method', diffusion_estimation_method, ...
        'msd_fit_input', msd_fit_input, ...
        'msd_fit_max_lag_s', msd_fit_max_lag_s, ...
        'msd_fit_max_track_fraction', msd_fit_max_track_fraction, ...
        'msd_fit_min_lag_span_s', msd_fit_min_lag_span_s, ...
        'msd_fit_min_r_squared', msd_fit_min_r_squared, ...
        'msd_fit_warn_r_squared', msd_fit_warn_r_squared, ...
        'msd_fit_Dt_method', "fixed condition-specific baseline (default zero)", ...
        'manual_exclusions', {manual_exclusions}, ...
        'repeat_exclusions', {repeat_exclusions}, ...
        'msd_max_lag_s', msd_max_lag_s, ...
        'msd_loglog', msd_loglog, ...
        'msd_sgolay_order', msd_sgolay_order, ...
        'msd_smoothing_window_frames', msd_smoothing_window_frames, ...
        'geometry_sgolay_order', geometry_sgolay_order, ...
        'geometry_smoothing_window_frames', geometry_smoothing_window_frames, ...
        'curvature_max_um_inv', curvature_max_um_inv, ...
        'curvature_plot_max_um_inv', curvature_plot_max_um_inv, ...
        'curvature_grid_points', curvature_grid_points, ...
        'curvature_bandwidth_um_inv', curvature_bandwidth, ...
        'dacf_max_lag_s', dacf_max_lag_s, ...
        'diffusion_translational_baseline_um2_s', diffusion_translational_baseline_um2_s, ...
        'diffusion_max_Dr_s_inv', diffusion_max_Dr_s_inv, ...
        'diffusion_min_fit_r_squared', diffusion_min_fit_r_squared, ...
        'diffusion_warn_fit_r_squared', diffusion_warn_fit_r_squared, ...
        'diffusion_grid_points', diffusion_grid_points, ...
        'diffusion_kde_method', "positive-support log-Gaussian KDE with 1/D Jacobian", ...
        'diffusion_bandwidth_log', diffusion_bandwidth_log, ...
        'diffusion_cv_bandwidth_limits_log', diffusion_kde_cv.search_limits_log, ...
        'diffusion_cv_grid_points', diffusion_cv_grid_points, ...
        'diffusion_kde_cv', diffusion_kde_cv, ...
        'diffusion_bandwidth_sensitivity_factors', diffusion_bandwidth_sensitivity_factors, ...
        'diffusion_bootstrap_samples', diffusion_bootstrap_samples, ...
        'diffusion_bootstrap_seed', diffusion_bootstrap_seed, ...
        'diffusion_plot_xlim_um2_s', diffusion_plot_xlim_um2_s, ...
        'diffusion_display_limits_um2_s', diffusion_display_limits, ...
        'diffusion_log_x', diffusion_log_x, ...
        'show_diffusion_peaks', show_diffusion_peaks, ...
        'show_diffusion_rug', show_diffusion_rug, ...
        'diffusion_omega_method', omega_method, ...
        'diffusion_fit_weights', fit_weights, ...
        'diffusion_pdf_band', "pointwise mean plus/minus 1 SEM across swimmer kernels; fixed bandwidth", ...
        'diffusion_kde_reference', "https://arxiv.org/abs/1804.08365", ...
        'diffusion_model_doi', "10.1103/PhysRevE.87.052305", ...
        'minimum_direction_step_um', minimum_direction_step_um);
    save(output_stem + "_results.mat", 'results', 'settings');
end

fprintf('\nSaved %s.png, %s.pdf, and %s.fig\n', ...
    output_stem, output_stem, output_stem);

%%%% LOCAL FUNCTIONS
function [trajectories, used_files, info] = load_curated_trajectories( ...
    data_dir, pixel_size_um, dt, minimum_duration_s, excluded_files)

    if ~isfolder(data_dir)
        error('Input folder not found: %s', data_dir);
    end

    files = dir(fullfile(data_dir, '*.csv'));
    files = files(~[files.isdir]);
    if isempty(files)
        error('No CSV files found in %s.', data_dir);
    end

    numeric_ids = nan(numel(files), 1);
    for k = 1:numel(files)
        [~, stem] = fileparts(files(k).name);
        numeric_ids(k) = str2double(stem);
    end
    [~, order] = sortrows([isnan(numeric_ids), numeric_ids], [1, 2]);
    files = files(order);

    trajectories = {};
    used_files = strings(0, 1);
    skipped_excluded = strings(0, 1);
    skipped_short = strings(0, 1);
    skipped_invalid = strings(0, 1);
    resampled_files = strings(0, 1);

    for k = 1:numel(files)
        file_name = string(files(k).name);
        if ismember(file_name, excluded_files)
            skipped_excluded(end+1,1) = file_name; %#ok<AGROW>
            continue;
        end

        file_path = fullfile(files(k).folder, files(k).name);
        T = readtable(file_path, 'VariableNamingRule', 'preserve');
        variables = string(T.Properties.VariableNames);

        if height(T) < 2 || ~all(ismember(["X", "Y"], variables))
            skipped_invalid(end+1,1) = file_name; %#ok<AGROW>
            continue;
        end

        x_px = double(T.X(:));
        y_px = double(T.Y(:));
        if ismember("FrameNumber", variables)
            frame_number = double(T.FrameNumber(:));
        else
            frame_number = (0:numel(x_px)-1)';
        end

        valid = isfinite(frame_number) & isfinite(x_px) & isfinite(y_px);
        frame_number = frame_number(valid);
        x_px = x_px(valid);
        y_px = y_px(valid);
        if numel(x_px) < 2
            skipped_invalid(end+1,1) = file_name; %#ok<AGROW>
            continue;
        end

        [frame_number, sort_order] = sort(frame_number);
        x_px = x_px(sort_order);
        y_px = y_px(sort_order);
        [frame_number, unique_index] = unique(frame_number, 'stable');
        x_px = x_px(unique_index);
        y_px = y_px(unique_index);

        t_native = (frame_number - frame_number(1)) * dt;
        x_um = (x_px - x_px(1)) * pixel_size_um;
        y_um = (y_px - y_px(1)) * pixel_size_um;

        if t_native(end) + 10*eps(t_native(end)) < minimum_duration_s
            skipped_short(end+1,1) = file_name; %#ok<AGROW>
            continue;
        end

        if any(abs(diff(t_native) - dt) > 100*eps(max(t_native(end), dt)))
            t_uniform = (0:floor(t_native(end) / dt))' * dt;
            x_um = interp1(t_native, x_um, t_uniform, 'pchip');
            y_um = interp1(t_native, y_um, t_uniform, 'pchip');
            resampled_files(end+1,1) = file_name; %#ok<AGROW>
        end

        if any(~isfinite(x_um)) || any(~isfinite(y_um))
            skipped_invalid(end+1,1) = file_name; %#ok<AGROW>
            continue;
        end

        trajectories{end+1,1} = [x_um(:), y_um(:)]; %#ok<AGROW>
        used_files(end+1,1) = file_name; %#ok<AGROW>
    end

    if isempty(trajectories)
        error('No valid trajectories remained in %s.', data_dir);
    end

    info = struct( ...
        'data_dir', string(data_dir), ...
        'excluded_files', skipped_excluded, ...
        'short_files', skipped_short, ...
        'invalid_files', skipped_invalid, ...
        'resampled_files', resampled_files);

    fprintf('Excluded repeats/outliers: %s\n', printable_list(skipped_excluded));
    fprintf('Skipped shorter than %.2f s: %s\n', ...
        minimum_duration_s, printable_list(skipped_short));
    fprintf('Skipped invalid/empty files: %s\n', printable_list(skipped_invalid));
    fprintf('Resampled files with frame gaps: %s\n', printable_list(resampled_files));
end

function output = printable_list(values)
    if isempty(values)
        output = 'none';
    else
        output = char(strjoin(values, ', '));
    end
end

function smoothed = smooth_trajectory_set(trajectories, order, window)
    if mod(window, 2) ~= 1 || window <= order
        error('Savitzky-Golay window must be odd and greater than its order.');
    end

    smoothed = cell(size(trajectories));
    for k = 1:numel(trajectories)
        xy = trajectories{k};
        if size(xy, 1) >= window
            x = sgolayfilt(xy(:,1), order, window);
            y = sgolayfilt(xy(:,2), order, window);
            smoothed{k} = [x, y];
        else
            smoothed{k} = xy;
        end
    end
end


function [tau, ensemble_msd, lower, upper, valid_track_count] = ...
    ensemble_geometric_msd(trajectories, dt, maximum_lag_s)

    maximum_lag = round(maximum_lag_s / dt);
    n_tracks = numel(trajectories);
    track_log_msd = nan(maximum_lag, n_tracks);
    pooled_log_sum = zeros(maximum_lag, 1);
    pooled_count = zeros(maximum_lag, 1);

    for k = 1:n_tracks
        xy = trajectories{k};
        n = size(xy, 1);
        for lag = 1:min(maximum_lag, n-1)
            displacement = xy(1+lag:n,:) - xy(1:n-lag,:);
            squared_displacement = sum(displacement.^2, 2);
            positive = squared_displacement(isfinite(squared_displacement) & ...
                squared_displacement > 0);
            if isempty(positive)
                continue;
            end
            log_values = log(positive);
            track_log_msd(lag,k) = mean(log_values);
            pooled_log_sum(lag) = pooled_log_sum(lag) + sum(log_values);
            pooled_count(lag) = pooled_count(lag) + numel(log_values);
        end
    end

    ensemble_msd = exp(pooled_log_sum ./ pooled_count);
    ensemble_msd(pooled_count == 0) = NaN;
    valid_track_count = sum(isfinite(track_log_msd), 2);
    log_sem = std(track_log_msd, 0, 2, 'omitnan') ./ ...
        sqrt(valid_track_count);
    log_sem(valid_track_count < 2) = NaN;
    lower = ensemble_msd .* exp(-log_sem);
    upper = ensemble_msd .* exp(log_sem);

    tau = (1:maximum_lag)' * dt;
    tau = [0; tau];
    ensemble_msd = [0; ensemble_msd];
    lower = [0; lower];
    upper = [0; upper];
    valid_track_count = [n_tracks; valid_track_count];
end


function curvature = menger_curvature(xy)
    n = size(xy, 1);
    curvature = nan(n, 1);
    if n < 3
        return;
    end

    p1 = xy(1:end-2,:);
    p2 = xy(2:end-1,:);
    p3 = xy(3:end,:);
    a = hypot(p2(:,1)-p1(:,1), p2(:,2)-p1(:,2));
    b = hypot(p3(:,1)-p2(:,1), p3(:,2)-p2(:,2));
    c = hypot(p3(:,1)-p1(:,1), p3(:,2)-p1(:,2));
    twice_area = abs((p2(:,1)-p1(:,1)).*(p3(:,2)-p1(:,2)) - ...
        (p2(:,2)-p1(:,2)).*(p3(:,1)-p1(:,1)));
    denominator = a .* b .* c;
    valid = isfinite(denominator) & denominator > eps;

    local_curvature = nan(size(denominator));
    local_curvature(valid) = 2 * twice_area(valid) ./ denominator(valid);
    curvature(2:end-1) = local_curvature;
end

function [pdf_mean, pdf_sem, pdf_matrix, valid_track_count] = ...
    ensemble_curvature_pdf(curvature_values, grid, bandwidth)

    n_tracks = numel(curvature_values);
    pdf_matrix = nan(n_tracks, numel(grid));
    for k = 1:n_tracks
        values = curvature_values{k};
        if numel(values) < 2
            continue;
        end
        reflected_values = [values(:); -values(:)];
        density = 2 * ksdensity(reflected_values, grid, ...
            'Bandwidth', bandwidth);
        area = trapz(grid, density);
        if isfinite(area) && area > 0
            pdf_matrix(k,:) = density / area;
        end
    end

    valid_tracks = all(isfinite(pdf_matrix), 2);
    valid_track_count = sum(valid_tracks);
    pdf_mean = mean(pdf_matrix(valid_tracks,:), 1, 'omitnan');
    pdf_sem = std(pdf_matrix(valid_tracks,:), 0, 1, 'omitnan') / ...
        sqrt(valid_track_count);
end

function peak_x = peak_location_parabolic(x, y)
    [~, index] = max(y);
    peak_x = x(index);
    if index <= 1 || index >= numel(x)
        return;
    end

    coefficients = polyfit(x(index-1:index+1), y(index-1:index+1), 2);
    candidate = -coefficients(2) / (2 * coefficients(1));
    if coefficients(1) < 0 && candidate >= x(index-1) && ...
            candidate <= x(index+1)
        peak_x = candidate;
    end
end

function [tau, ensemble_mean, ensemble_sem, correlation_matrix, pair_counts] = ...
    ensemble_directional_autocorrelation(trajectories, dt, ...
    maximum_lag_s, minimum_step_um)

    maximum_lag = round(maximum_lag_s / dt);
    n_tracks = numel(trajectories);
    correlation_matrix = nan(maximum_lag + 1, n_tracks);
    pair_counts = zeros(maximum_lag + 1, n_tracks);

    for k = 1:n_tracks
        displacement = diff(trajectories{k}, 1, 1);
        step_length = hypot(displacement(:,1), displacement(:,2));
        direction = nan(size(displacement));
        valid_step = isfinite(step_length) & step_length > minimum_step_um;
        direction(valid_step,:) = displacement(valid_step,:) ./ ...
            step_length(valid_step);

        n = size(direction, 1);
        for lag = 0:min(maximum_lag, n-1)
            first = direction(1:n-lag,:);
            second = direction(1+lag:n,:);
            valid_pair = all(isfinite(first), 2) & all(isfinite(second), 2);
            pair_counts(lag+1,k) = sum(valid_pair);
            if any(valid_pair)
                products = sum(first(valid_pair,:) .* second(valid_pair,:), 2);
                correlation_matrix(lag+1,k) = mean(products);
            end
        end
    end

    tau = (0:maximum_lag)' * dt;
    ensemble_mean = mean(correlation_matrix, 2, 'omitnan');
    valid_count = sum(isfinite(correlation_matrix), 2);
    ensemble_sem = std(correlation_matrix, 0, 2, 'omitnan') ./ ...
        sqrt(valid_count);
    ensemble_sem(valid_count < 2) = NaN;
end

function fit_result = fit_damped_dacf(tau, correlation, sem)
    valid = tau > 0 & isfinite(tau) & isfinite(correlation);
    fit_tau = tau(valid);
    fit_correlation = correlation(valid);
    fit_sem = sem(valid);

    positive_sem = fit_sem(isfinite(fit_sem) & fit_sem > 0);
    if isempty(positive_sem)
        sem_floor = 0.05;
    else
        sem_floor = max(0.02, 0.25 * median(positive_sem));
    end
    weights = 1 ./ max(fit_sem, sem_floor).^2;
    weights(~isfinite(weights)) = 1;
    weights = weights / mean(weights);

    zero_index = find(fit_correlation <= 0, 1, 'first');
    if isempty(zero_index)
        omega_initial = 2*pi / max(2*fit_tau(end), eps);
    else
        first_zero_time = fit_tau(zero_index);
        omega_initial = pi / max(2*first_zero_time, eps);
    end
    Dr_initial = 0.5;

    model = @(q, t) exp(-exp(q(1))*t) .* cos(exp(q(2))*t);
    objective = @(q) sum(weights .* ...
        (fit_correlation - model(q, fit_tau)).^2);
    options = optimset('Display', 'off', 'MaxFunEvals', 4000, ...
        'MaxIter', 2000, 'TolX', 1e-10, 'TolFun', 1e-10);
    fitted_q = fminsearch(objective, log([Dr_initial, omega_initial]), ...
        options);

    Dr = exp(fitted_q(1));
    omega = exp(fitted_q(2));
    fitted_curve = exp(-Dr*tau) .* cos(omega*tau);
    fitted_values = model(fitted_q, fit_tau);
    residual_sum_squares = sum((fit_correlation - fitted_values).^2);
    total_sum_squares = sum((fit_correlation - mean(fit_correlation)).^2);
    r_squared = 1 - residual_sum_squares / total_sum_squares;

    fit_result = struct( ...
        'Dr_s_inv', Dr, ...
        'decay_time_s', 1/Dr, ...
        'omega_rad_s', omega, ...
        'period_s', 2*pi/omega, ...
        'r_squared', r_squared, ...
        'fitted_curve', fitted_curve);
end

function [parameters, fits] = swimmer_msd_fit_parameters( ...
    trajectories, files, dt, maximum_lag_s, maximum_track_fraction, ...
    minimum_lag_span_s, Dt, maximum_Dr, minimum_r_squared, ...
    warning_r_squared)

    if maximum_track_fraction <= 0 || maximum_track_fraction > 1 || ...
            maximum_lag_s <= 0 || minimum_lag_span_s <= 0
        error('Invalid single-swimmer MSD fit window settings.');
    end
    n_tracks = numel(trajectories);
    parameters = table(string(files(:)), 'VariableNames', {'file'});
    numeric_fields = ["duration_s", "fit_max_lag_s", "fit_lag_count", ...
        "fit_speed_um_s", "omega_rad_s", "Dr_s_inv", "Dt_um2_s", ...
        "Deff_um2_s", "fit_r_squared", "fit_rmse_um2"];
    for name = numeric_fields
        parameters.(name) = nan(n_tracks, 1);
    end
    logical_fields = ["included_in_diffusion_pdf", "weak_fit", ...
        "noise_decay_observed", "rotation_observed", ...
        "omega_at_upper_bound", "Dr_at_lower_bound", ...
        "Dr_at_upper_bound"];
    for name = logical_fields
        parameters.(name) = false(n_tracks, 1);
    end
    parameters.fit_status = repmat("not fitted", n_tracks, 1);
    parameters.exclusion_reason = repmat("", n_tracks, 1);
    fits = cell(n_tracks, 1);

    for k = 1:n_tracks
        xy = trajectories{k};
        n = size(xy,1);
        maximum_lag = min(round(maximum_lag_s/dt), ...
            floor(maximum_track_fraction*(n-1)));
        parameters.duration_s(k) = (n-1)*dt;
        parameters.fit_max_lag_s(k) = maximum_lag*dt;
        parameters.Dt_um2_s(k) = Dt;
        if maximum_lag*dt < minimum_lag_span_s || maximum_lag < 20
            parameters.fit_status(k) = "fit window too short";
            parameters.exclusion_reason(k) = parameters.fit_status(k);
            continue;
        end

        % Algebraic, not geometric: for each lag average |r(t+tau)-r(t)|^2
        % over ALL available time origins of this ONE particle.
        [tau, observed, pair_count] = single_particle_algebraic_msd( ...
            xy, dt, maximum_lag);
        valid = tau > 0 & isfinite(observed) & pair_count > 0;
        tau = tau(valid);
        observed = observed(valid);
        pair_count = pair_count(valid);
        parameters.fit_lag_count(k) = numel(tau);
        if numel(tau) < 20
            parameters.fit_status(k) = "too few MSD lags";
            parameters.exclusion_reason(k) = parameters.fit_status(k);
            continue;
        end

        fit = fit_circle_swimmer_msd(tau, observed, Dt, dt, maximum_Dr);
        fits{k} = struct('lag_s', tau, 'observed_um2', observed, ...
            'model_um2', fit.fitted_curve, 'origin_count', pair_count);
        parameters.fit_speed_um_s(k) = fit.v_um_s;
        parameters.omega_rad_s(k) = fit.omega_rad_s;
        parameters.Dr_s_inv(k) = fit.Dr_s_inv;
        parameters.fit_r_squared(k) = fit.r_squared;
        parameters.fit_rmse_um2(k) = fit.rmse_um2;
        parameters.fit_status(k) = fit.status;
        parameters.weak_fit(k) = fit.r_squared < warning_r_squared;
        parameters.Dr_at_lower_bound(k) = fit.Dr_s_inv < 1e-5;
        parameters.Dr_at_upper_bound(k) = fit.Dr_s_inv >= maximum_Dr*(1-1e-5);
        parameters.omega_at_upper_bound(k) = ...
            fit.omega_rad_s >= (pi/dt)*(1-1e-5);
        parameters.noise_decay_observed(k) = ...
            fit.Dr_s_inv*max(tau) >= 1;
        parameters.rotation_observed(k) = ...
            fit.omega_rad_s*max(tau) >= 2*pi;
        Deff = effective_circle_diffusion( ...
            fit.v_um_s, fit.omega_rad_s, fit.Dr_s_inv, Dt);
        parameters.Deff_um2_s(k) = Deff;

        if fit.status ~= "ok"
            parameters.exclusion_reason(k) = fit.status;
        elseif parameters.Dr_at_upper_bound(k) || ...
                parameters.omega_at_upper_bound(k)
            parameters.exclusion_reason(k) = "MSD fit reached search bound";
        elseif ~isfinite(Deff) || Deff <= 0
            parameters.exclusion_reason(k) = "no positive finite long-time diffusion";
        elseif ~isfinite(fit.r_squared) || fit.r_squared < minimum_r_squared
            parameters.exclusion_reason(k) = "MSD incompatible with circle-swimmer model";
        else
            parameters.included_in_diffusion_pdf(k) = true;
        end
    end
end

function [tau, msd, origin_count] = single_particle_algebraic_msd(xy, dt, maximum_lag)
    n = size(xy,1);
    maximum_lag = min(maximum_lag, n-1);
    tau = (1:maximum_lag)'*dt;
    msd = nan(maximum_lag,1);
    origin_count = zeros(maximum_lag,1);
    for lag = 1:maximum_lag
        displacement = xy(1+lag:n,:) - xy(1:n-lag,:);
        squared = sum(displacement.^2,2);
        valid = isfinite(squared);
        origin_count(lag) = sum(valid);
        if origin_count(lag) > 0
            msd(lag) = mean(squared(valid));
        end
    end
end

function fit = fit_circle_swimmer_msd(tau, observed, Dt, dt, maximum_Dr)
    % Profile out v^2 analytically at every (omega, Dr): the remaining
    % least-squares search is only 2-D, with nonnegative parameters.
    % A coarse scan plus multistart bounded refinement avoids dependence on
    % a single initial period or decay rate.
    omega_bound = pi/dt;
    omegas = unique(min(omega_bound, ...
        [0, 0.25, 0.5, 0.8, 1.3, 2, 3.2, 5, 8, 13, 21, 35]));
    Drs = unique(min(maximum_Dr, ...
        [0, 0.015, 0.05, 0.15, 0.4, 1, 2.5, 6, 15, 40]));
    scale = max(max(observed), 1);
    scores = inf(numel(omegas), numel(Drs));
    for i = 1:numel(omegas)
        for j = 1:numel(Drs)
            predicted = circle_msd_profile( ...
                tau, observed, omegas(i), Drs(j), Dt);
            scores(i,j) = sum(((predicted-observed)/scale).^2);
        end
    end
    [best_by_omega, best_Dr_index] = min(scores, [], 2);
    [~, omega_order] = sort(best_by_omega);
    selected = omega_order(1:min(7,numel(omega_order)));
    options = optimoptions('lsqnonlin', 'Display', 'off', ...
        'MaxIterations', 180, 'MaxFunctionEvaluations', 700, ...
        'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10);
    best_score = Inf;
    best_q = [NaN, NaN];
    best_exit_flag = -1;
    for i = selected(:)'
        seed = [omegas(i), Drs(best_Dr_index(i))];
        objective = @(q) (circle_msd_profile( ...
            tau, observed, q(1), q(2), Dt)-observed)/scale;
        [q, ~, residual, exit_flag] = lsqnonlin( ...
            objective, seed, [0,0], [omega_bound,maximum_Dr], options);
        score = sum(residual.^2);
        if score < best_score
            best_score = score;
            best_q = q;
            best_exit_flag = exit_flag;
        end
    end
    [predicted, speed_squared] = circle_msd_profile( ...
        tau, observed, best_q(1), best_q(2), Dt);
    residual = predicted-observed;
    fit = struct('v_um_s', sqrt(speed_squared), ...
        'omega_rad_s', best_q(1), 'Dr_s_inv', best_q(2), ...
        'r_squared', fit_r_squared(sum(residual.^2), ...
            sum((observed-mean(observed)).^2)), ...
        'rmse_um2', sqrt(mean(residual.^2)), ...
        'status', "ok", 'fitted_curve', predicted);
    if best_exit_flag <= 0
        fit.status = "MSD least-squares fit did not converge";
    end
end

function [predicted, speed_squared] = circle_msd_profile( ...
    tau, observed, omega, Dr, Dt)
    active_unit = circle_swimmer_msd(tau, 1, omega, Dr, 0);
    active_target = observed-4*Dt*tau;
    denominator = sum(active_unit.^2);
    if denominator <= realmin
        speed_squared = 0;
    else
        speed_squared = max(0, sum(active_unit.*active_target)/denominator);
    end
    predicted = 4*Dt*tau+speed_squared*active_unit;
end

function msd = circle_swimmer_msd(tau, v, omega, Dr, Dt)
    % Algebraically identical to Marine et al. Eq. (5). The series branch
    % keeps the Dr=omega=0 ballistic limit and very short lags accurate.
    tau = tau(:);
    z = Dr-1i*omega;
    zt = z*tau;
    integrated_correlation = zeros(size(tau));
    small = abs(zt) < 0.01;
    t = tau(small);
    integrated_correlation(small) = t.^2 .* ...
        (1/2 - z*t/6 + z.^2.*t.^2/24 - z.^3.*t.^3/120 + ...
        z.^4.*t.^4/720);
    if any(~small)
        integrated_correlation(~small) = ...
            (exp(-zt(~small))-1+zt(~small))/z^2;
    end
    msd = 4*Dt*tau + 2*v^2*real(integrated_correlation);
end

function [parameters, fitted_curves] = swimmer_diffusion_parameters( ...
    trajectories, files, curvature_values, dt, minimum_step_um, tau, ...
    correlation_matrix, pair_counts, Dt, maximum_Dr, minimum_r_squared, ...
    warning_r_squared)

    n_tracks = numel(trajectories);
    parameters = table(string(files(:)), 'VariableNames', {'file'});
    numeric_fields = ["duration_s", "mean_speed_um_s", ...
        "mean_curvature_um_inv", "mean_signed_curvature_um_inv", ...
        "mean_signed_omega_rad_s", "omega_rad_s", "Dr_s_inv", ...
        "Dt_um2_s", "Deff_um2_s", "fit_r_squared", ...
        "fit_weighted_r_squared", "fit_rmse", "valid_step_fraction", ...
        "valid_turn_count", "fit_lag_count"];
    for name = numeric_fields
        parameters.(name) = nan(n_tracks, 1);
    end
    logical_fields = ["included_in_diffusion_pdf", "weak_fit", ...
        "noise_decay_observed", "Dr_at_lower_bound", "Dr_at_upper_bound"];
    for name = logical_fields
        parameters.(name) = false(n_tracks, 1);
    end
    parameters.fit_status = repmat("not fitted", n_tracks, 1);
    parameters.exclusion_reason = repmat("", n_tracks, 1);
    fitted_curves = nan(size(correlation_matrix));

    for k = 1:n_tracks
        xy = trajectories{k};
        displacement = diff(xy, 1, 1);
        step_length = hypot(displacement(:,1), displacement(:,2));
        valid_speed = isfinite(step_length);
        valid_step = valid_speed & step_length > minimum_step_um;
        valid_turn = valid_step(1:end-1) & valid_step(2:end);

        parameters.duration_s(k) = (size(xy,1)-1) * dt;
        parameters.mean_speed_um_s(k) = mean(step_length(valid_speed)) / dt;
        parameters.mean_curvature_um_inv(k) = mean(curvature_values{k});
        parameters.Dt_um2_s(k) = Dt;
        parameters.valid_step_fraction(k) = mean(valid_step);
        parameters.valid_turn_count(k) = sum(valid_turn);
        if sum(valid_turn) < 3
            parameters.exclusion_reason(k) = "too few valid heading increments";
            continue;
        end

        first = displacement(1:end-1,:);
        second = displacement(2:end,:);
        % Wrapped signed increments avoid an unwrap jump across a masked
        % near-stationary step. Do NOT connect directions across such gaps.
        turns = atan2(first(:,1).*second(:,2)-first(:,2).*second(:,1), ...
            sum(first .* second, 2));
        turns = turns(valid_turn);
        omega_signed = mean(turns) / dt;
        omega = abs(omega_signed); % chirality sign does not affect D_eff
        local_arc_length = (step_length(1:end-1) + step_length(2:end)) / 2;
        parameters.mean_signed_curvature_um_inv(k) = ...
            sum(turns) / sum(local_arc_length(valid_turn));
        parameters.mean_signed_omega_rad_s(k) = omega_signed;
        parameters.omega_rad_s(k) = omega;

        % Each fit uses ONLY this swimmer's DACF. The ensemble fit in panel
        % l is never used for a swimmer's diffusion coefficient. Longer-lag
        % estimates with fewer valid time origins receive less weight.
        fit = fit_swimmer_Dr(tau, correlation_matrix(:,k), ...
            pair_counts(:,k), omega, maximum_Dr);
        parameters.Dr_s_inv(k) = fit.Dr_s_inv;
        parameters.fit_r_squared(k) = fit.r_squared;
        parameters.fit_weighted_r_squared(k) = fit.weighted_r_squared;
        parameters.fit_rmse(k) = fit.rmse;
        parameters.fit_lag_count(k) = fit.lag_count;
        parameters.fit_status(k) = fit.status;
        parameters.Dr_at_lower_bound(k) = fit.at_lower_bound;
        parameters.Dr_at_upper_bound(k) = fit.at_upper_bound;
        fitted_curves(:,k) = fit.fitted_curve;
        parameters.weak_fit(k) = fit.r_squared < warning_r_squared;
        parameters.noise_decay_observed(k) = ...
            fit.Dr_s_inv * parameters.duration_s(k) >= 1;
        Deff = effective_circle_diffusion(parameters.mean_speed_um_s(k), ...
            omega, fit.Dr_s_inv, Dt);
        parameters.Deff_um2_s(k) = Deff;

        if fit.status ~= "ok"
            parameters.exclusion_reason(k) = fit.status;
        elseif fit.at_upper_bound
            parameters.exclusion_reason(k) = "Dr search bound reached";
        elseif ~isfinite(Deff) || Deff < 0
            parameters.exclusion_reason(k) = "no finite diffusive long-time limit";
        elseif ~isfinite(fit.r_squared) || fit.r_squared < minimum_r_squared
            parameters.exclusion_reason(k) = "DACF incompatible with fixed mean-turn model";
        else
            parameters.included_in_diffusion_pdf(k) = true;
        end
    end
end

function fit = fit_swimmer_Dr(tau, correlation, pair_count, omega, maximum_Dr)
    fit = struct('Dr_s_inv', NaN, 'r_squared', NaN, ...
        'weighted_r_squared', NaN, 'rmse', NaN, 'lag_count', 0, ...
        'at_lower_bound', false, 'at_upper_bound', false, ...
        'status', "insufficient DACF data", 'fitted_curve', nan(size(tau)));
    valid = tau > 0 & isfinite(tau) & isfinite(correlation) & ...
        isfinite(pair_count) & pair_count > 0;
    fit.lag_count = sum(valid);
    if fit.lag_count < 4 || ~isfinite(omega) || omega < 0 || ...
            ~isfinite(maximum_Dr) || maximum_Dr <= 0
        return;
    end
    t = tau(valid);
    y = correlation(valid);
    weights = pair_count(valid) / mean(pair_count(valid));
    model = @(Dr) exp(-Dr*t) .* cos(omega*t);
    objective = @(Dr) sum(weights .* (y-model(Dr)).^2);

    % A coarse global scan followed by local refinement avoids dependence
    % on an initial Dr and resolves the noiseless Dr=0 boundary explicitly.
    search_grid = unique([0, linspace(0, maximum_Dr, 100), ...
        logspace(-6, log10(maximum_Dr), 120)]);
    objectives = arrayfun(objective, search_grid);
    [~, best] = min(objectives);
    lower = search_grid(max(1, best-1));
    upper = search_grid(min(numel(search_grid), best+1));
    options = optimset('Display', 'off', 'TolX', 1e-10);
    [candidate, ~, exit_flag] = fminbnd(objective, lower, upper, options);
    if exit_flag <= 0
        fit.status = "Dr minimization did not converge";
        return;
    end
    candidates = [0, search_grid(best), candidate, maximum_Dr];
    [~, best] = min(arrayfun(objective, candidates));
    Dr = candidates(best);
    residual = y-model(Dr);
    total_ss = sum((y-mean(y)).^2);
    weighted_mean = sum(weights.*y) / sum(weights);
    weighted_total_ss = sum(weights.*(y-weighted_mean).^2);
    fit.Dr_s_inv = Dr;
    fit.r_squared = fit_r_squared(sum(residual.^2), total_ss);
    fit.weighted_r_squared = fit_r_squared( ...
        sum(weights.*residual.^2), weighted_total_ss);
    fit.rmse = sqrt(mean(residual.^2));
    fit.at_lower_bound = Dr == 0;
    fit.at_upper_bound = maximum_Dr-Dr <= 1e-6*maximum_Dr;
    fit.status = "ok";
    fit.fitted_curve = exp(-Dr*tau) .* cos(omega*tau);
end

function r_squared = fit_r_squared(residual_ss, total_ss)
    if total_ss <= eps
        if residual_ss <= 100*eps
            r_squared = 1;
        else
            r_squared = -Inf;
        end
    else
        r_squared = 1-residual_ss/total_ss;
    end
end

function Deff = effective_circle_diffusion(v, omega, Dr, Dt)
    if any(~isfinite([v, omega, Dr, Dt])) || any([v, omega, Dr, Dt] < 0)
        Deff = NaN;
    elseif v == 0
        Deff = Dt;
    elseif Dr == 0 && omega == 0
        % A noiseless straight swimmer remains ballistic, not diffusive.
        Deff = Inf;
    else
        Deff = Dt + v^2*Dr / (2*(Dr^2 + omega^2));
    end
end

function [pdf_mean, pdf_sem, kernels] = swimmer_diffusion_pdf(values, grid, bandwidth)
    % One log-normal kernel per swimmer, normalized on 0<D<infinity.
    % The 1/D Jacobian is essential: this is p(D), not p(log(D)).
    % Kernels remain zero at D=0; no reflection or epsilon pseudodata.
    % Actual zero coefficients would require a separately modeled atom.
    values = values(:);
    grid = grid(:)';
    if isempty(values) || any(~isfinite(values)) || any(values <= 0)
        error('analysis_2jklm:NonpositiveDiffusion', ...
            ['Log-KDE requires strictly positive coefficients. Genuine zeros ', ...
            'need a separate probability mass, not an added small constant.']);
    end
    if ~isscalar(bandwidth) || ~isfinite(bandwidth) || bandwidth <= 0 || ...
            any(~isfinite(grid))
        error('Log-KDE bandwidth must be positive and the grid must be finite.');
    end
    kernels = zeros(numel(values), numel(grid));
    positive = grid > 0;
    standardized = (log(grid(positive))-log(values))/bandwidth;
    kernels(:,positive) = exp(-0.5*standardized.^2) ./ ...
        (sqrt(2*pi)*bandwidth*grid(positive));
    pdf_mean = mean(kernels, 1);
    pdf_sem = std(kernels, 0, 1) / sqrt(numel(values));
    % Retain the kernel SEM as a numerical diagnostic. Panel m instead
    % plots bootstrap percentile bands, conditional on h and model inputs.
end

function [bandwidth, cv] = select_log_kde_bandwidth(groups, search_limits, n_grid)
    n_conditions = numel(groups);
    all_values = vertcat(groups{:});
    for c = 1:n_conditions
        values = groups{c};
        if numel(values) < 2 || any(~isfinite(values)) || any(values <= 0)
            error('Each log-KDE group needs at least two positive coefficients.');
        end
    end
    log_values = log(all_values);
    if range(log_values) <= 100*eps(max(1, max(abs(log_values))))
        error('A density bandwidth cannot be identified from identical coefficients.');
    end
    [~, ~, reference] = ksdensity(log_values);
    if isempty(search_limits)
        search_limits = [max(1e-4, 0.05*reference), max(1e-3, 4*reference)];
    end
    if numel(search_limits) ~= 2 || any(~isfinite(search_limits)) || ...
            search_limits(1) <= 0 || search_limits(2) <= search_limits(1) || ...
            ~isscalar(n_grid) || n_grid < 5 || n_grid ~= round(n_grid)
        error('Invalid log-bandwidth search bounds or grid size.');
    end
    log_candidates = linspace(log(search_limits(1)), log(search_limits(2)), n_grid);
    candidates = exp(log_candidates);
    scores = nan(1, n_grid);
    condition_scores = nan(n_conditions, n_grid);
    for j = 1:n_grid
        [scores(j), condition_scores(:,j)] = log_kde_cv_score(groups, candidates(j));
    end
    [~, best] = max(scores);
    lower = log_candidates(max(1, best-1));
    upper = log_candidates(min(n_grid, best+1));
    options = optimset('Display', 'off', 'TolX', 1e-9);
    refined = fminbnd(@(q) -log_kde_cv_score(groups, exp(q)), lower, upper, options);
    choices = [search_limits(1), candidates(best), exp(refined), search_limits(2)];
    choice_scores = arrayfun(@(h) log_kde_cv_score(groups, h), choices);
    [best_score, best] = max(choice_scores);
    bandwidth = choices(best);
    at_bound = min(abs(log(bandwidth)-log(search_limits))) < 1e-6;
    if at_bound
        warning('analysis_2jklm:BandwidthSearchBound', ...
            'Selected KDE bandwidth reaches a search bound; inspect the CV diagnostic.');
    end
    cv = struct('bandwidth_log', candidates, 'balanced_log_score', scores, ...
        'condition_log_scores', condition_scores, 'selected_bandwidth_log', bandwidth, ...
        'selected_log_score', best_score, 'search_limits_log', search_limits, ...
        'normal_reference_bandwidth_log', reference, 'at_search_bound', at_bound, ...
        'criterion', "mean within-condition held-out log p(D), equal condition weights");
end

function [score, condition_scores] = log_kde_cv_score(groups, bandwidth)
    condition_scores = nan(numel(groups), 1);
    for c = 1:numel(groups)
        z = log(groups{c}(:));
        n = numel(z);
        terms = -0.5*((z-z')/bandwidth).^2;
        terms(1:n+1:end) = -Inf; % leave this swimmer out completely
        row_max = max(terms, [], 2);
        log_sum = row_max+log(sum(exp(terms-row_max), 2));
        held_out_log_density = log_sum-log(n-1)-log(bandwidth) - ...
            0.5*log(2*pi)-z; % -z = log of the 1/D Jacobian
        condition_scores(c) = mean(held_out_log_density);
    end
    % Equal condition weighting prevents the slightly larger group from
    % setting the smoothing. Neighbors always come from the SAME group.
    score = mean(condition_scores);
end

function [lower, upper, peak_samples] = bootstrap_kde(kernels, grid, n_bootstrap, seed)
    if ~isscalar(n_bootstrap) || n_bootstrap < 2 || n_bootstrap ~= round(n_bootstrap)
        error('Use an integer of at least two bootstrap replicates.');
    end
    % Resample independent swimmers, not trajectory frames. Restore the
    % caller's random stream so these diagnostics cannot alter analyses.
    state = rng;
    restore_rng = onCleanup(@() rng(state));
    rng(seed, 'twister');
    n = size(kernels, 1);
    indices = randi(n, n, n_bootstrap);
    columns = repelem((1:n_bootstrap)', n);
    weights = accumarray([indices(:), columns], 1, [n, n_bootstrap])/n;
    bootstrap_density = weights'*kernels;
    percentiles = prctile(bootstrap_density, [2.5, 97.5], 1);
    lower = percentiles(1,:);
    upper = percentiles(2,:);
    peak_samples = nan(n_bootstrap, 1);
    for b = 1:n_bootstrap
        peak_samples(b) = peak_location_parabolic(grid, bootstrap_density(b,:));
    end
    % These are conditional sampling ranges for the smoothed estimate.
    % They do not propagate fitted Dr/v/omega errors, nor bandwidth choice;
    % the latter is checked separately by the sensitivity curves.
end

function fig = plot_single_swimmer_msd_fit_diagnostics(results, colors, names)
    fig = figure('Color', 'w', 'Position', [80, 80, 1350, 750]);
    layout = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', ...
        'Padding', 'loose');
    sgtitle(layout, ['Single-swimmer algebraic MSD: observed (points) ', ...
        'and Eq. (5) fit (line)'], 'FontSize', 15);
    quality_levels = [10, 50, 90];
    for c = 1:numel(results)
        parameters = results(c).swimmer_parameters;
        included = find(parameters.included_in_diffusion_pdf);
        quality = parameters.fit_r_squared(included);
        for column = 1:numel(quality_levels)
            target = prctile(quality, quality_levels(column));
            [~, order] = min(abs(quality-target));
            k = included(order);
            fit = results(c).diffusion_msd_fits{k};
            ax = nexttile(layout, (c-1)*numel(quality_levels)+column);
            hold(ax, 'on');
            plot(ax, fit.lag_s, fit.observed_um2, '.', ...
                'Color', colors(c,:), 'MarkerSize', 8, ...
                'DisplayName', 'Observed');
            plot(ax, fit.lag_s, fit.model_um2, '-', ...
                'Color', 0.6*colors(c,:), 'LineWidth', 2, ...
                'DisplayName', 'Analytical fit');
            xlim(ax, [0, max(fit.lag_s)]);
            ylim(ax, padded_zero_ylim( ...
                [fit.observed_um2; fit.model_um2], 0.08));
            title(ax, sprintf('%s | %s | R^2 = %.2f', ...
                names(c), parameters.file(k), ...
                parameters.fit_r_squared(k)), 'FontWeight', 'normal');
            text(ax, 0.97, 0.93, sprintf('D_{eff} = %.1f \\mum^2/s', ...
                parameters.Deff_um2_s(k)), 'Units', 'normalized', ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
                'FontSize', 10, 'Color', [0.3, 0.3, 0.3]);
            if ~parameters.noise_decay_observed(k)
                text(ax, 0.97, 0.83, 'decay not resolved', ...
                    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                    'VerticalAlignment', 'top', 'FontSize', 10, ...
                    'Color', [0.65, 0.2, 0.2]);
            end
            xlabel(ax, '\tau (s)', 'Interpreter', 'tex');
            ylabel(ax, 'MSD (\mum^2)', 'Interpreter', 'tex');
            set(ax, 'FontName', 'Arial', 'FontSize', 11, ...
                'LineWidth', 1.1, 'Box', 'off', 'TickDir', 'out');
            if c == 1 && column == 1
                legend(ax, 'Location', 'northwest', 'Box', 'off');
            end
        end
    end
end

function fig = plot_log_kde_diagnostics(results, cv, colors, names, factors, limits)
    fig = figure('Color', 'w', 'Position', [80, 80, 1350, 900]);
    layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'loose');
    styles = ["--", "-", ":"];
    for c = 1:numel(results)
        ax = nexttile(layout, c); hold(ax, 'on');
        for f = 1:numel(factors)
            plot(ax, results(c).diffusion_grid_um2_s, ...
                results(c).diffusion_sensitivity_pdf(f,:), ...
                'LineStyle', styles(mod(f-1, numel(styles))+1), ...
                'Color', colors(c,:), 'LineWidth', 1.8, ...
                'DisplayName', sprintf('h = %.3f (%.2fx)', ...
                cv.used_bandwidth_log*factors(f), factors(f)));
        end
        xlim(ax, limits);
        ylim(ax, padded_zero_ylim(results(c).diffusion_sensitivity_pdf(:), 0.1));
        v = results(c).diffusion_values_um2_s;
        v = v(v >= limits(1) & v <= limits(2))';
        plot(ax, [v; v], [zeros(size(v)); 0.025*ax.YLim(2)*ones(size(v))], ...
            'Color', colors(c,:), 'HandleVisibility', 'off');
        xlabel(ax, 'D_{eff} (\mum^2/s)', 'Interpreter', 'tex');
        ylabel(ax, 'PDF (s/\mum^2)', 'Interpreter', 'tex');
        title(ax, names(c)+" (bandwidth sensitivity)", 'FontWeight', 'normal');
        legend(ax, 'Location', 'northeast', 'Box', 'off', 'FontSize', 10);
        set(ax, 'FontName', 'Arial', 'FontSize', 12, 'LineWidth', 1.2, ...
            'TickDir', 'out', 'Box', 'off');
    end

    ax = nexttile(layout, 3); hold(ax, 'on');
    near = cv.bandwidth_log >= cv.used_bandwidth_log/3 & ...
        cv.bandwidth_log <= cv.used_bandwidth_log*3;
    if ~any(near)
        near(:) = true;
    end
    for c = 1:numel(results)
        plot(ax, cv.bandwidth_log(near), cv.condition_log_scores(c,near), ...
            'Color', colors(c,:), 'LineWidth', 1.6, 'DisplayName', names(c));
    end
    plot(ax, cv.bandwidth_log(near), cv.balanced_log_score(near), ...
        'k-', 'LineWidth', 2, 'DisplayName', 'Balanced score');
    xline(ax, cv.used_bandwidth_log, 'k--', 'HandleVisibility', 'off');
    set(ax, 'XScale', 'log', 'FontName', 'Arial', 'FontSize', 12, ...
        'LineWidth', 1.2, 'TickDir', 'out', 'Box', 'off');
    xlabel(ax, 'Log-space bandwidth h');
    ylabel(ax, 'Mean held-out log density');
    title(ax, 'Leave-one-swimmer-out bandwidth selection', 'FontWeight', 'normal');
    legend(ax, 'Location', 'southeast', 'Box', 'off', 'FontSize', 10);

    ax = nexttile(layout, 4); hold(ax, 'on');
    all_peaks = vertcat(results.diffusion_bootstrap_peak_um2_s);
    maximum = max(150, 1.05*max(all_peaks));
    edges = linspace(0, maximum, 41);
    for c = 1:numel(results)
        histogram(ax, results(c).diffusion_bootstrap_peak_um2_s, edges, ...
            'Normalization', 'pdf', 'FaceColor', colors(c,:), ...
            'FaceAlpha', 0.30, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        interval = results(c).diffusion_peak_ci95_um2_s;
        text(ax, 0.97, 0.92-0.09*(c-1), sprintf('%s: %.1f--%.1f', ...
            names(c), interval), 'Units', 'normalized', ...
            'HorizontalAlignment', 'right', 'Color', colors(c,:), 'FontSize', 10);
    end
    xlim(ax, [0, maximum]);
    set(ax, 'FontName', 'Arial', 'FontSize', 12, 'LineWidth', 1.2, ...
        'TickDir', 'out', 'Box', 'off');
    xlabel(ax, 'Bootstrap PDF peak (\mum^2/s)', 'Interpreter', 'tex');
    ylabel(ax, 'PDF of bootstrap peaks');
    title(ax, 'Peak variability (selected bandwidth held fixed)', 'FontWeight', 'normal');
end

function test_circle_msd_fit()
    tau = (0:150)'/15;
    v = 12;
    omega = 1.4;
    Dr = 0.35;
    analytic = circle_swimmer_msd(tau, v, omega, Dr, 0);
    assert(analytic(1) == 0 && all(analytic >= -1e-9), ...
        'Circle-swimmer MSD must start at zero and remain nonnegative.');
    ballistic = circle_swimmer_msd(tau, v, 0, 0, 0.3);
    assert(max(abs(ballistic-(v^2*tau.^2+1.2*tau))) < 1e-8, ...
        'The ballistic-plus-translational limit is incorrect.');
    pure_circle = circle_swimmer_msd(tau, v, omega, 0, 0);
    assert(max(abs(pure_circle-2*v^2*(1-cos(omega*tau))/omega^2)) < 1e-8, ...
        'The noiseless circular limit is incorrect.');
    long_tau = 1e5;
    long_msd = circle_swimmer_msd(long_tau, v, omega, Dr, 0.3);
    assert(abs(long_msd/(4*long_tau)- ...
        effective_circle_diffusion(v,omega,Dr,0.3)) < 0.02, ...
        'The long-time MSD slope must agree with the paper diffusivity.');
    straight_xy = [v*tau, zeros(size(tau))];
    [lags, measured, origins] = single_particle_algebraic_msd( ...
        straight_xy, 1/15, 100);
    assert(max(abs(measured-v^2*lags.^2)) < 1e-8 && ...
        isequal(origins, (numel(tau)-(1:100))'), ...
        'Single-particle MSD must average all squared-displacement origins.');
    fit = fit_circle_swimmer_msd(tau(2:end), analytic(2:end), 0, 1/15, 75);
    assert(fit.status == "ok" && abs(fit.v_um_s-v) < 0.01 && ...
        abs(fit.omega_rad_s-omega) < 0.001 && ...
        abs(fit.Dr_s_inv-Dr) < 0.001 && fit.r_squared > 1-1e-8, ...
        'The single-swimmer least-squares fit failed synthetic recovery.');
    fprintf('Algebraic MSD and analytical circle-swimmer fit self-tests passed.\n');
end

function test_diffusion_model()
    assert(abs(effective_circle_diffusion(10, 0, 2, 0.2)-25.2) < 1e-12, ...
        'The noncircling 2-D active-Brownian limit must be v^2/(2*Dr).');
    assert(effective_circle_diffusion(10, 2, 0, 0.2) == 0.2, ...
        'A noiseless circle adds no long-time active diffusion.');
    assert(isinf(effective_circle_diffusion(10, 0, 0, 0)), ...
        'A ballistic swimmer must not be assigned a finite diffusion.');
    assert(effective_circle_diffusion(0, 0, 0, 0.2) == 0.2);

    tau = (0:60)'/15;
    counts = 300-(0:60)';
    for Dr = [0, 0.02, 0.3, 1, 4]
        C = exp(-Dr*tau).*cos(1.7*tau);
        fit = fit_swimmer_Dr(tau, C, counts, 1.7, 75);
        assert(abs(fit.Dr_s_inv-Dr) < 1e-6 && fit.r_squared > 1-1e-10, ...
            'Individual DACF fits failed to recover a known decay.');
    end
    grid = [0, logspace(-4, 3, 3000)];
    [density, sem, kernels] = swimmer_diffusion_pdf([1; 2; 4], grid, 0.4);
    assert(all(kernels >= 0, 'all') && all(sem >= 0));
    assert(abs(trapz(grid, density)-1) < 1e-4, ...
        'The log-KDE must normalize after applying the 1/D Jacobian.');
    assert(all(kernels(:,1) == 0), 'Log-normal kernels must vanish at D=0.');
    assert(size(kernels,1) == 3, 'There must be one kernel per swimmer.');
    single_density = swimmer_diffusion_pdf(10, grid, 0.4);
    assert(abs(peak_location_parabolic(grid, single_density)-10*exp(-0.4^2)) < 0.01, ...
        'A single kernel must recover the analytic log-normal mode.');
    assert(abs(trapz(grid, grid.*single_density)-10*exp(0.4^2/2)) < 1e-3, ...
        'A single kernel must recover the analytic log-normal mean.');
    rejected_zero = false;
    try
        swimmer_diffusion_pdf([0; 1], grid, 0.4);
    catch exception
        rejected_zero = strcmp(exception.identifier, 'analysis_2jklm:NonpositiveDiffusion');
    end
    assert(rejected_zero, 'A true zero must not be replaced with epsilon pseudodata.');
    groups = {exp([-0.7; -0.4; -0.1; 0.2; 0.8]); exp([0.1; 0.3; 0.6; 0.9; 1.2])};
    h = select_log_kde_bandwidth(groups, [0.03, 1.5], 41);
    h_scaled = select_log_kde_bandwidth(cellfun(@(v) v*100, groups, ...
        'UniformOutput', false), [0.03, 1.5], 41);
    assert(abs(h-h_scaled) < 1e-5, 'CV bandwidth must be invariant to diffusion units.');
    [lo, hi, peaks] = bootstrap_kde(kernels, grid, 40, 123);
    [lo_again, hi_again, peaks_again] = bootstrap_kde(kernels, grid, 40, 123);
    assert(isequal(lo, lo_again) && isequal(hi, hi_again) && isequal(peaks, peaks_again));
    assert(all(lo <= hi) && all(peaks > 0));
    fprintf('Diffusion model, DACF recovery, log-KDE, CV, and bootstrap self-tests passed.\n');
end

function plot_band(ax, x, lower, upper, color, alpha)
    x = x(:);
    lower = lower(:);
    upper = upper(:);
    valid = isfinite(x) & isfinite(lower) & isfinite(upper);
    if ~any(valid)
        return;
    end

    x = x(valid);
    lower = lower(valid);
    upper = upper(valid);
    fill(ax, [x; flipud(x)], [lower; flipud(upper)], color, ...
        'FaceAlpha', alpha, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
end

function plot_msd_panel(ax, results, colors, names, use_loglog, dt, max_lag_s)
    for c = 1:numel(results)
        x = results(c).msd_time_s;
        y = results(c).msd_um2 / 1e3;
        lower = results(c).msd_lower_um2 / 1e3;
        upper = results(c).msd_upper_um2 / 1e3;

        if use_loglog
            plot_mask = x > 0 & y > 0 & lower > 0 & upper > 0 & ...
                isfinite(x) & isfinite(y) & isfinite(lower) & isfinite(upper);
        else
            plot_mask = isfinite(x) & isfinite(y) & ...
                isfinite(lower) & isfinite(upper);
        end

        plot_band(ax, x(plot_mask), lower(plot_mask), upper(plot_mask), ...
            colors(c,:), 0.18);
        plot(ax, x(plot_mask), y(plot_mask), ...
            'Color', colors(c,:), 'LineWidth', 2.2, ...
            'DisplayName', names(c));
    end

    if use_loglog
        set(ax, 'XScale', 'log', 'YScale', 'log');
        xlim(ax, [dt, max_lag_s]);
    else
        xlim(ax, [0, max_lag_s]);
        ylim(ax, padded_zero_ylim([results.msd_upper_um2] / 1e3, 0.05));
    end

    xlabel(ax, '\tau (s)', 'Interpreter', 'tex');
    ylabel(ax, 'MSD (10^3 \mum^2)', 'Interpreter', 'tex');
    legend(ax, 'Location', 'northwest', 'Box', 'off');
end

function plot_optional_dacf_fit(ax, x, fitted_curve, color, show_fit)
    if ~show_fit
        return;
    end
    plot(ax, x, fitted_curve, '--', 'Color', color, 'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
end

function limits = padded_zero_ylim(values, padding_fraction)
    values = values(isfinite(values));
    if isempty(values)
        limits = [0, 1];
        return;
    end
    maximum = max(values);
    limits = [0, maximum * (1 + padding_fraction)];
end

function style_panel(ax, letter)
    set(ax, 'FontName', 'Arial', 'FontSize', 12, 'LineWidth', 1.2, ...
        'TickDir', 'out', 'Box', 'off', 'Layer', 'top');
    ax.XColor = [0, 0, 0];
    ax.YColor = [0, 0, 0];
    text(ax, -0.12, 1.02, letter, 'Units', 'normalized', ...
        'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
        'Clipping', 'off');
end
