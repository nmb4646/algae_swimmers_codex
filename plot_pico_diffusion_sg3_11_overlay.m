clear; clc; close all;

% One-panel comparison of the SG3/11 diffusion PDFs. Both methods were
% computed after identical trajectory smoothing in the abc comparison.
source_mat = "compare_pico_diffusion_methods_abc_sg1_3_results.mat";
output_stem = "pico_diffusion_sg3_11_bc_overlay";
display_xlim_um2_s = [0, 40]; % Peak region; PDFs normalize over full support.
colors = [0.12, 0.47, 0.71; 1.00, 0.25, 0.39]; % 0%, 2% mucin
condition_names = ["0% mucin"; "2% mucin"];

assert(isfile(source_mat), ...
    'Run compare_pico_diffusion_methods_abc.m to produce the sweep results.');
S = load(source_mat, 'settings', 'group_index', 'b_values', 'c_values', ...
    'x_grid');
variant = find(S.settings.sg_order == 3 & ...
    S.settings.sg_window_frames == 11);
assert(isscalar(variant), 'SG3/11 is absent from the saved sweep.');
h = S.settings.kde_bandwidth_log;
x = S.x_grid(:)';
visible = x >= display_xlim_um2_s(1) & ...
    x <= display_xlim_um2_s(2);

fig = figure('Color', 'w', 'Position', [80, 80, 1100, 720]);
ax = axes(fig);
hold(ax, 'on');
peak_um2_s = nan(2, 2); % rows: 0%, 2%; columns: b, c
maximum_visible_pdf = 0;

for condition = 1:2
    idx = S.group_index{condition};
    coefficient_sets = {S.b_values(variant, idx), ...
        S.c_values(variant, idx)};
    for method = 1:2
        density = positive_log_kde(coefficient_sets{method}, x, h);
        [~, peak_idx] = max(density);
        peak_um2_s(condition, method) = x(peak_idx);
        maximum_visible_pdf = max(maximum_visible_pdf, ...
            max(density(visible)));
        if method == 1
            style = '--';
            method_name = 'b: DACF';
        else
            style = '-';
            method_name = 'c: MSD fit';
        end
        plot(ax, x(visible), density(visible), style, ...
            'Color', colors(condition,:), 'LineWidth', 2.8, ...
            'DisplayName', sprintf('%s, %s', ...
            condition_names(condition), method_name));
    end
end

xlim(ax, display_xlim_um2_s);
ylim(ax, [0, 1.10*maximum_visible_pdf]);
xlabel(ax, 'Modeled D_{eff} (\mum^2/s)', 'Interpreter', 'tex');
ylabel(ax, 'Probability density (s/\mum^2)', 'Interpreter', 'tex');
title(ax, 'Long-time diffusion after SG3/11 trajectory smoothing', ...
    'FontWeight', 'normal', 'FontSize', 17);
legend(ax, 'Location', 'northeast', 'Box', 'off', 'FontSize', 12);
set(ax, 'FontName', 'Arial', 'FontSize', 14, 'LineWidth', 1.5, ...
    'TickDir', 'out', 'Box', 'off');
ax.Toolbar.Visible = 'off';
drawnow;

exportgraphics(fig, output_stem + ".png", 'Resolution', 400);
exportgraphics(fig, output_stem + ".pdf", 'ContentType', 'vector');
savefig(fig, output_stem + ".fig");

fprintf('SG3/11 b/c modes (um^2/s):\n');
for condition = 1:2
    fprintf('%s: %.3f / %.3f (n = %d)\n', ...
        condition_names(condition), peak_um2_s(condition,:), ...
        numel(S.group_index{condition}));
end
fprintf('Saved %s.png and %s.pdf\n', output_stem, output_stem);

function density = positive_log_kde(values, x, bandwidth)
    density = zeros(size(x));
    positive = x > 0;
    z = (log(x(positive)) - log(values(:))) / bandwidth;
    density(positive) = mean(exp(-0.5*z.^2) / ...
        (sqrt(2*pi)*bandwidth), 1) ./ x(positive);
end
