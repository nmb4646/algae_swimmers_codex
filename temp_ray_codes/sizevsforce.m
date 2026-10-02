% =========================================================================
% Drag Force vs. Particle Diameter via Mesh Activation Model (chi = 500 nm)
% X-axis range: 0.2 - 1.6 um; Fitted/simulated lines plotted as dashed lines ('--')
% =========================================================================

clear; clc; close all;

% 1. Experimental Input Data from Tables
d_data = [1.0, 1.2, 1.5];                       % Particle diameters (um)

% 0% Mucin Data: [Mean Speed, Std Dev] (um/s)
v0_data  = [57.820, 55.189, 48.314];
std_v0   = [4.089,  3.671,  5.135];

% 2% Mucin Data: [Mean Speed, Std Dev] (um/s)
v2_data  = [52.635, 44.179, 32.000];
std_v2   = [3.271,  4.506,  4.480];

% 2. Physical Constants & Parameters
eta0 = 1.0e-3;                                  % Water viscosity (Pa·s)
chi  = 0.5;                                     % Mesh size chi = 500 nm = 0.5 um

% 3. Friction Coefficient Factors (pN / (um/s))
f0_coef = 3 * pi * 1e-3 * d_data;               % Pure water friction factor
f2_coef = f0_coef .* exp(d_data ./ chi);        % Mesh friction factor

% 4. Drag Force & Error Bar Calculations (pN)
F0_data = f0_coef .* v0_data;                   % 0% Mucin Force
F2_data = f2_coef .* v2_data;                   % 2% Mucin Force

err_F0  = f0_coef .* std_v0;                    % Propagated Force Error (0%)
err_F2  = f2_coef .* std_v2;                    % Propagated Force Error (2%)

% 5. Smooth Simulated Curve Fitting (X-axis range: 0.2 - 1.6 um)
d_fine = linspace(0.2, 1.6, 200);

p0 = polyfit(d_data, F0_data, 1);
F0_fit = polyval(p0, d_fine);

p2 = polyfit(d_data, F2_data, 2);
F2_fit = polyval(p2, d_fine);

% 6. Plotting Initialization
figure('Color', 'w', 'Position', [200, 200, 480, 400]);
hold on;

% Academic Palette Matching Reference Image
color_blue = [0.02, 0.45, 0.75];                % Deep Blue (0% Mucin)
color_pink = [0.95, 0.20, 0.40];                % Vibrant Pink-Red (2% Mucin)

% --- 0% Mucin Data Points & Dashed Simulated Line ---
h_e0 = errorbar(d_data, F0_data, err_F0, 'o', ...
    'Color', color_blue, 'MarkerFaceColor', color_blue, ...
    'MarkerSize', 8, 'LineWidth', 2.0, 'CapSize', 6, ...
    'DisplayName', '0% Mucin Data');

h_f0 = plot(d_fine, F0_fit, '--', 'Color', color_blue, ...
    'LineWidth', 3.0, 'DisplayName', '0% Mucin Simulated Curve');

% --- 2% Mucin Data Points & Dashed Simulated Line ---
h_e2 = errorbar(d_data, F2_data, err_F2, 's', ...
    'Color', color_pink, 'MarkerFaceColor', color_pink, ...
    'MarkerSize', 8, 'LineWidth', 2.0, 'CapSize', 6, ...
    'DisplayName', '2% Mucin (\chi = 500 nm)');

h_f2 = plot(d_fine, F2_fit, '--', 'Color', color_pink, ...
    'LineWidth', 3.0, 'DisplayName', '2% Mucin Simulated Curve');

% 7. Minimalist Academic Styling (X-axis range: 0.2 - 1.6 um)
box off; grid off;                               % Clean canvas

set(gca, 'LineWidth', 1.8, ...                  % Thick main axis lines
         'TickDir', 'out', ...                  % Ticks pointing outwards
         'FontName', 'Arial', ...
         'FontSize', 12, ...
         'FontWeight', 'bold', ...
         'XColor', 'k', 'YColor', 'k');

xlabel('Particle Diameter (\mum)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Microscopic Drag Force (pN)', 'FontSize', 13, 'FontWeight', 'bold');

xlim([0.2, 1.6]);                               % X-axis set to 0.2 - 1.6 um
ylim([0, 11]);

% Borderless Legend
legend([h_e0, h_f0, h_e2, h_f2], 'Location', 'northwest', ...
    'FontSize', 10, 'Box', 'off');