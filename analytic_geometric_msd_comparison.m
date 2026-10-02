clear; clc; close all;

%%%% ANALYTIC GEOMETRIC-MSD COMPARISON FOR THE PICOALGAE
% The usual circle-swimmer MSD is E[R^2]. The geometric MSD used in Figure
% 2j is exp(E[log(R^2)]). This script evaluates the latter without random
% trajectory simulations, using the circle swimmer's exact displacement
% characteristic function and a one-dimensional Bessel/log integral.
%
% The parameters below are a *diagnostic* parameterization, not a fitted
% long-time forecast: v is the arithmetic mean per-swimmer speed, omega and
% Dr are from the existing ensemble DACF fit, and Dt is an adjustable
% translational-noise baseline. A finite Dt also regularizes the logarithm
% near returns to the same position.
%
% Theory: Kurzthaler & Franosch, Soft Matter 13, 6396 (2017), DOI
% 10.1039/c7sm00873b, derive the circle swimmer characteristic function.
% The logarithmic Bessel identity below converts it to a geometric MSD.

%%%% SETTINGS
input_mat = 'analysis_2jklm_results.mat';
output_stem = 'analytic_geometric_msd_Dr_comparison';
translational_diffusion_um2_s = 0.214; % illustrative 1-um-radius sphere in water
paper_sphere_radius_um = 1.0;     % matches D2.m; test actual cell size separately
paper_temperature_K = 293;
paper_viscosity_pa_s = [0.001; 0.001]; % 2% entry is a WATER counterfactual
angular_fourier_modes = 50;           % retain modes -50:50
wave_number_points = 300;
wave_number_limits_um_inv = [1e-5, 45];
maximum_lag_s = 25;
image_dpi = 300;

if ~isfile(input_mat)
    error('Run analysis_2jklm.m first to create %s.', input_mat);
end
loaded = load(input_mat, 'results', 'settings');
results = loaded.results;
settings = loaded.settings;
dt = settings.dt_s;
maximum_lag = min(round(maximum_lag_s/dt), ...
    numel(results(1).msd_time_s)-1);
lag_s = (1:maximum_lag)' * dt;
n_conditions = numel(results);
if numel(paper_viscosity_pa_s) ~= n_conditions
    error('Provide one viscosity value for each condition.');
end
colors = [0.12, 0.47, 0.71; 1.00, 0.25, 0.39];

comparison = repmat(struct(), n_conditions, 1);
for c = 1:n_conditions
    data_dir = fullfile(settings.base_dir, settings.condition_dirs(c));
    files = string(results(c).files);
    [raw_tracks, speeds] = load_original_tracks( ...
        data_dir, files, settings.pixel_size_um, dt, ...
        settings.msd_sgolay_order, settings.msd_smoothing_window_frames);

    v = mean(speeds);
    omega = results(c).dacf_fit.omega_rad_s;
    Dr = results(c).dacf_fit.Dr_s_inv;
    Dt = translational_diffusion_um2_s;
    [model_geometric, model_arithmetic] = circle_geometric_msd( ...
        lag_s, v, omega, Dr, Dt, angular_fourier_modes, ...
        wave_number_points, wave_number_limits_um_inv);
    boltzmann_J_K = 1.380649e-23;
    a_m = paper_sphere_radius_um * 1e-6;
    eta = paper_viscosity_pa_s(c);
    paper_Dr = boltzmann_J_K * paper_temperature_K / ...
        (8*pi*eta*a_m^3);
    paper_Dt = boltzmann_J_K * paper_temperature_K / ...
        (6*pi*eta*a_m) * 1e12;
    [paper_geometric, paper_arithmetic] = circle_geometric_msd( ...
        lag_s, v, omega, paper_Dr, paper_Dt, ...
        angular_fourier_modes, wave_number_points, ...
        wave_number_limits_um_inv);

    comparison(c).condition = results(c).condition;
    comparison(c).lag_s = lag_s;
    comparison(c).model_geometric_um2 = model_geometric;
    comparison(c).model_arithmetic_um2 = model_arithmetic;
    comparison(c).paper_geometric_um2 = paper_geometric;
    comparison(c).paper_arithmetic_um2 = paper_arithmetic;
    comparison(c).experimental_smoothed_geometric_um2 = ...
        results(c).msd_um2(2:maximum_lag+1);
    comparison(c).experimental_raw_geometric_um2 = ...
        pooled_geometric_msd(raw_tracks, maximum_lag);
    comparison(c).v_um_s = v;
    comparison(c).omega_rad_s = omega;
    comparison(c).Dr_s_inv = Dr;
    comparison(c).Dt_um2_s = Dt;
    comparison(c).paper_Dr_s_inv = paper_Dr;
    comparison(c).paper_Dt_um2_s = paper_Dt;
    comparison(c).paper_sphere_radius_um = paper_sphere_radius_um;
    comparison(c).paper_viscosity_pa_s = eta;
    comparison(c).n_swimmers = numel(raw_tracks);

    last_index = maximum_lag;
    fprintf(['%s: n=%d, v=%.3f um/s, omega=%.3f rad/s, ', ...
        'Dr=%.3f s^-1, Dt=%.3f um^2/s\n'], ...
        string(comparison(c).condition), comparison(c).n_swimmers, ...
        v, omega, Dr, Dt);
    fprintf(['  Paper-style thermal Dr=%.4f s^-1; ', ...
        'Dt=%.4f um^2/s (radius %.2f um, eta %.3g Pa s)\n'], ...
        paper_Dr, paper_Dt, paper_sphere_radius_um, eta);
    fprintf(['  At %.1f s: DACF analytic=%.1f, thermal analytic=%.1f, ', ...
        'raw experimental=%.1f, smoothed experimental=%.1f um^2\n'], ...
        lag_s(last_index), ...
        model_geometric(last_index), ...
        paper_geometric(last_index), ...
        comparison(c).experimental_raw_geometric_um2(last_index), ...
        comparison(c).experimental_smoothed_geometric_um2(last_index));
end

%%%% ANALYTIC LIMIT CHECKS
% In two dimensions, E[R^2]=4Dt and geometric MSD=e^(-EulerGamma)*4Dt.
test_lag_s = [0.2; 1; 5];
[brownian_geometric, brownian_arithmetic] = circle_geometric_msd( ...
    test_lag_s, 0, 0, 1, translational_diffusion_um2_s, 1, ...
    60, wave_number_limits_um_inv);
relative_test_error = max(abs( ...
    brownian_geometric ./ brownian_arithmetic - exp(-euler_gamma)));
assert(relative_test_error < 1e-8, ...
    'Analytic geometric-MSD evaluator failed its Brownian limit test.');
% A perfect circle plus translational Brownian noise has a noncentral
% chi-square radial law. Its geometric MSD is m^2 exp(E1(m^2/(4Dt*t))),
% where m=2(v/omega)sin(omega*t/2). This checks the Fourier/log integral,
% rather than only the arithmetic-MSD expression.
test_circle_time_s = [0.5; 1.0];
test_circle_v_um_s = 45;
test_circle_omega_rad_s = 1.4;
test_circle_m2 = 4*(test_circle_v_um_s/test_circle_omega_rad_s)^2 ...
    .* sin(test_circle_omega_rad_s*test_circle_time_s/2).^2;
expected_circle_geometric = test_circle_m2 .* exp(expint( ...
    test_circle_m2 ./ (4*translational_diffusion_um2_s*test_circle_time_s)));
[computed_circle_geometric, ~] = circle_geometric_msd( ...
    test_circle_time_s, test_circle_v_um_s, ...
    test_circle_omega_rad_s, 0, translational_diffusion_um2_s, ...
    angular_fourier_modes, wave_number_points, ...
    wave_number_limits_um_inv);
relative_circle_test_error = max(abs(computed_circle_geometric ./ ...
    expected_circle_geometric - 1));
assert(relative_circle_test_error < 0.01, ...
    'Analytic geometric-MSD evaluator failed its perfect-circle limit test.');

%%%% COMPARISON FIGURE
fig = figure('Color', 'w', 'Position', [80, 100, 1120, 470]);
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle(layout, sprintf(['Thermal-sphere comparison: a = %.1f \\mum, ', ...
    '\\eta = %.1f mPa s in both conditions (2%% counterfactual)'], ...
    paper_sphere_radius_um, 1e3*paper_viscosity_pa_s(1)), ...
    'FontWeight', 'normal', 'FontSize', 11);
for c = 1:n_conditions
    ax = nexttile(layout, c);
    hold(ax, 'on');
    plot(ax, lag_s, comparison(c).experimental_raw_geometric_um2, ...
        ':', 'Color', colors(c,:), 'LineWidth', 1.8, ...
        'DisplayName', 'Experiment, raw positions');
    plot(ax, lag_s, comparison(c).experimental_smoothed_geometric_um2, ...
        '-', 'Color', colors(c,:), 'LineWidth', 2.1, ...
        'DisplayName', 'Experiment, Fig. 2j preprocessing');
    plot(ax, lag_s, comparison(c).model_geometric_um2, ...
        '-.', 'Color', [0.42, 0.42, 0.42], 'LineWidth', 1.8, ...
        'DisplayName', 'Model: DACF decay as D_r');
    plot(ax, lag_s, comparison(c).paper_geometric_um2, ...
        '--', 'Color', [0.08, 0.08, 0.08], 'LineWidth', 2.1, ...
        'DisplayName', 'Model: thermal sphere D_r');
    set(ax, 'XScale', 'log', 'YScale', 'log', 'FontName', 'Arial', ...
        'FontSize', 12, 'LineWidth', 1.2, 'Box', 'on');
    xlabel(ax, '\tau (s)', 'Interpreter', 'tex');
    if c == 1
        ylabel(ax, 'Geometric MSD (\mum^2)', 'Interpreter', 'tex');
    end
    title(ax, string(comparison(c).condition), 'FontWeight', 'normal');
    xlim(ax, [lag_s(1), lag_s(end)]);
    grid(ax, 'on');
    legend(ax, 'Location', 'northwest', 'Box', 'off', 'FontSize', 9);
end
exportgraphics(fig, output_stem + ".png", 'Resolution', image_dpi);
exportgraphics(fig, output_stem + ".pdf", 'ContentType', 'vector');
save(output_stem + ".mat", 'comparison', 'relative_test_error', ...
    'relative_circle_test_error', ...
    'paper_sphere_radius_um', 'paper_temperature_K', ...
    'paper_viscosity_pa_s', ...
    'angular_fourier_modes', 'wave_number_points', ...
    'wave_number_limits_um_inv');

%%%% LOCAL FUNCTIONS
function [tracks, speeds] = load_original_tracks(data_dir, files, ...
    pixel_size_um, dt, sg_order, sg_window)
    tracks = cell(numel(files), 1);
    speeds = nan(numel(files), 1);
    for j = 1:numel(files)
        T = readtable(fullfile(data_dir, files(j)), ...
            'VariableNamingRule', 'preserve');
        frame = double(T.FrameNumber(:));
        x = double(T.X(:));
        y = double(T.Y(:));
        valid = isfinite(frame) & isfinite(x) & isfinite(y);
        frame = frame(valid); x = x(valid); y = y(valid);
        [frame, order] = sort(frame);
        x = x(order); y = y(order);
        [frame, unique_index] = unique(frame, 'stable');
        x = x(unique_index); y = y(unique_index);
        t = (frame - frame(1)) * dt;
        x = (x - x(1)) * pixel_size_um;
        y = (y - y(1)) * pixel_size_um;
        if any(abs(diff(t) - dt) > 100*eps(max(t(end), dt)))
            uniform_t = (0:floor(t(end)/dt))' * dt;
            x = interp1(t, x, uniform_t, 'pchip');
            y = interp1(t, y, uniform_t, 'pchip');
        end
        tracks{j} = [x(:), y(:)];
        smooth_x = sgolayfilt(x, sg_order, sg_window);
        smooth_y = sgolayfilt(y, sg_order, sg_window);
        speeds(j) = mean(hypot(gradient(smooth_x)/dt, ...
            gradient(smooth_y)/dt));
    end
end

function geometric = pooled_geometric_msd(tracks, maximum_lag)
    log_sum = zeros(maximum_lag, 1);
    count = zeros(maximum_lag, 1);
    for j = 1:numel(tracks)
        xy = tracks{j};
        n = size(xy, 1);
        for lag = 1:min(maximum_lag, n-1)
            displacement = xy(lag+1:n,:) - xy(1:n-lag,:);
            squared = sum(displacement.^2, 2);
            squared = squared(isfinite(squared) & squared > 0);
            log_sum(lag) = log_sum(lag) + sum(log(squared));
            count(lag) = count(lag) + numel(squared);
        end
    end
    geometric = exp(log_sum ./ count);
    geometric(count == 0) = NaN;
end

function [geometric, arithmetic] = circle_geometric_msd(t, v, omega, ...
    Dr, Dt, mode_count, k_points, k_bounds)
    % The Langevin model is dr=v(cos theta,sin theta)dt+sqrt(2Dt)dW,
    % dtheta=omega dt+sqrt(2Dr)dB. Its characteristic function F(k,t)
    % equals E[J0(kR)] after averaging over initial headings.
    %
    % For an isotropic 2-D Gaussian with the same arithmetic MSD M(t),
    % F_ref(k,t)=exp(-M(t)k^2/4) and E[log R_ref^2]=log M-gamma.
    % The Bessel identity then yields the *exact* geometric MSD:
    % log G=log M-gamma+2 integral_0^infty (F_ref-F) dk/k.
    % Here F is evaluated by Fourier modes of the exact angular generator.
    t = t(:);
    if any(t <= 0) || Dt <= 0
        error('All lags and translational diffusivity must be positive.');
    end
    lambda = -Dr + 1i*omega;
    if abs(lambda) < 1e-12
        arithmetic = 4*Dt*t + v^2*t.^2;
    else
        arithmetic = 4*Dt*t + 2*v^2*real( ...
            (exp(lambda*t) - 1 - lambda*t) / lambda^2);
    end
    if v == 0
        geometric = exp(-euler_gamma) * arithmetic;
        return;
    end

    k = logspace(log10(k_bounds(1)), log10(k_bounds(2)), k_points)';
    modes = (-mode_count:mode_count)';
    center = mode_count + 1;
    n_modes = numel(modes);
    initial = zeros(n_modes, 1);
    initial(center) = 1;
    characteristic = zeros(numel(k), numel(t));

    diagonal = -Dr*modes.^2 + 1i*omega*modes;
    for j = 1:numel(k)
        off_diagonal = 1i*k(j)*v/2 * ones(n_modes-1, 1);
        generator = diag(diagonal) + diag(off_diagonal, 1) + ...
            diag(off_diagonal, -1);
        [eigenvectors, eigenvalues] = eig(generator, 'vector');
        coefficients = eigenvectors \ initial;
        weights = eigenvectors(center,:).' .* coefficients;
        active_characteristic = real(weights.' * ...
            exp(eigenvalues * t.'));
        characteristic(j,:) = active_characteristic .* ...
            exp(-Dt*k(j)^2*t.');
    end

    gaussian_reference = exp(-k.^2 * arithmetic.' / 4);
    correction = 2 * trapz(log(k), ...
        gaussian_reference - characteristic, 1).';
    geometric = exp(log(arithmetic) - euler_gamma + correction);
end

function value = euler_gamma()
    value = 0.5772156649015328606;
end
