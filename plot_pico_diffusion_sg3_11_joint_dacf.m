clear; clc; close all;

% Exploratory SG3/11 PDF: fit each swimmer's DACF with both D_r and omega
% free. Current fixed-omega and method-c PDFs are shown for context.
base_dir = string(fileparts(mfilename('fullpath')));
S = load(fullfile(base_dir, ...
    "compare_pico_diffusion_methods_abc_sg1_3_results.mat"));
P = readtable(fullfile(base_dir, ...
    "analysis_2jklm_swimmer_parameters.csv"), 'TextType', 'string');
B = load(fullfile(base_dir, "analysis_2jklm_results.mat"), 'results');
variant = find(S.settings.sg_order == 3 & ...
    S.settings.sg_window_frames == 11);
assert(isscalar(variant));
dt = 1/15;
h = S.settings.kde_bandwidth_log;
x = S.x_grid(S.x_grid > 0);
folders = ["picoalgae 0% mucin"; "picoalgae 2% mucin"];
method_names = ["current"; "max 2 s"; "skip first 4 lags"; ...
    "free amplitude"; "free omega"];
joint_values = cell(2,1);

for condition = 1:2
    cohort = S.group_index{condition};
    count = numel(cohort);
    Drs = nan(count, 5);
    omegas = nan(count, 5);
    amplitudes = nan(count, 5);
    r2s = nan(count, 5);
    speeds = S.speed_um_s(variant, cohort)';
    measured_D = S.b_values(variant, cohort)';

    for swimmer = 1:count
        index = cohort(swimmer);
        path = fullfile(base_dir, "data/Picoalgae in 0-2% mucin", ...
            folders(condition), P.file(index));
        T = readtable(path, 'VariableNamingRule', 'preserve');
        frames = double(T.FrameNumber(:));
        assert(all(diff(frames) == 1), 'Unexpected frame gap.');
        xy = [double(T.X(:)), double(T.Y(:))];
        xy = (xy - xy(1,:))*0.48;
        xy = [sgolayfilt(xy(:,1),3,11), ...
            sgolayfilt(xy(:,2),3,11)];
        [t, y, weights, omega] = track_dacf(xy,dt);
        weights = weights/mean(weights);
        baseline = S.dacf_Dr_s_inv(variant,index);
        [check, amplitude] = fit_dr(t,y,weights,omega,false);
        assert(abs(check-baseline) < 1e-6, ...
            'Diagnostic does not reproduce the saved Dr.');
        Drs(swimmer,1) = baseline;
        omegas(swimmer,1) = omega;
        amplitudes(swimmer,1) = amplitude;
        r2s(swimmer,1) = fit_r2(t,y,baseline,omega,1);

        [Drs(swimmer,2), amplitudes(swimmer,2)] = ...
            fit_dr(t(t<=2),y(t<=2),weights(t<=2),omega,false);
        omegas(swimmer,2) = omega;
        r2s(swimmer,2) = fit_r2(t,y,Drs(swimmer,2),omega,1);

        later = t >= 5*dt;
        [Drs(swimmer,3), amplitudes(swimmer,3)] = ...
            fit_dr(t(later),y(later),weights(later),omega,false);
        omegas(swimmer,3) = omega;
        r2s(swimmer,3) = fit_r2(t,y,Drs(swimmer,3),omega,1);

        [Drs(swimmer,4), amplitudes(swimmer,4)] = ...
            fit_dr(t,y,weights,omega,true);
        omegas(swimmer,4) = omega;
        r2s(swimmer,4) = fit_r2(t,y,Drs(swimmer,4), ...
            omega,amplitudes(swimmer,4));

        ensemble_omega = B.results(condition).dacf_fit.omega_rad_s;
        [Drs(swimmer,5), omegas(swimmer,5)] = ...
            fit_dr_free_omega(t,y,weights,baseline,omega,ensemble_omega);
        amplitudes(swimmer,5) = 1;
        r2s(swimmer,5) = fit_r2(t,y,Drs(swimmer,5), ...
            omegas(swimmer,5),1);
    end

    fprintf('\nCondition %d%%, n=%d\n', (condition-1)*2, count);
    for method = 1:5
        D = speeds.^2 .* Drs(:,method) ./ ...
            (2*(Drs(:,method).^2 + omegas(:,method).^2));
        density = mean(exp(-0.5*((log(x)-log(D))/h).^2) / ...
            (sqrt(2*pi)*h), 1) ./ x;
        [~, peak] = max(density);
        fprintf(['%s: Dr median %.3f IQR [%.3f %.3f]; ', ...
            'R2 median %.3f; amp median %.3f; omega median %.3f; ', ...
            'D mode %.3f; D median %.3f; Dr ratio to current %.3f\n'], ...
            method_names(method), median(Drs(:,method)), ...
            prctile(Drs(:,method),[25,75]), median(r2s(:,method)), ...
            median(amplitudes(:,method)),median(omegas(:,method)), ...
            x(peak),median(D),median(Drs(:,method)./Drs(:,1)));
    end
    fprintf('Original saved method-b D max absolute discrepancy %.3g\n', ...
        max(abs(measured_D - speeds.^2 .* Drs(:,1) ./ ...
        (2*(Drs(:,1).^2 + omegas(:,1).^2)))));
    delta_omega = abs(omegas(:,5)-omegas(:,1));
    poor = r2s(:,1)<0.5;
    fprintf(['Free-omega median |delta omega| %.3f rad/s, ', ...
        'median 4-s phase shift %.3f rad, ', ...
        'median R2 improvement %.3f\n'], ...
        median(delta_omega), 4*median(delta_omega), ...
        median(r2s(:,5)-r2s(:,1)));
    fprintf(['Poor current fits: n=%d, median Dr current/free ', ...
        '%.3f/%.3f, median R2 current/free %.3f/%.3f\n'], ...
        sum(poor), median(Drs(poor,1)), median(Drs(poor,5)), ...
        median(r2s(poor,1)), median(r2s(poor,5)));
    durations = P.duration_s(cohort);
    fprintf('Duration median poor/good: %.2f/%.2f s; Spearman duration vs current R2 %.3f\n', ...
        median(durations(poor)), median(durations(~poor)), ...
        corr(durations,r2s(:,1),'Type','Spearman'));
    joint_values{condition} = speeds.^2 .* Drs(:,5) ./ (2*(Drs(:,5).^2 + omegas(:,5).^2));
end

% Same positive log-KDE as the saved a/b/c comparison; only the DACF fit
% changes. Display the peak region without renormalizing after cropping.
condition_names = ["0% mucin";"2% mucin"];
colors = [0.12,0.47,0.71;1.00,0.25,0.39];
fig = figure('Color','w','Position',[80,80,1300,590]);
layout = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
visible = x<=40;
for condition = 1:2
    cohort = S.group_index{condition};
    values = {S.b_values(variant,cohort),joint_values{condition}, ...
        S.c_values(variant,cohort)};
    curves = cell(1,3);
    modes = nan(1,3);
    for method = 1:3
        sample = values{method}(:);
        z = (log(x)-log(sample))/h;
        curves{method} = mean(exp(-0.5*z.^2)/(sqrt(2*pi)*h),1)./x;
        [~,peak_index] = max(curves{method});
        modes(method) = x(peak_index);
    end
    ax = nexttile(layout,condition);
    hold(ax,'on');
    plot(ax,x(visible),curves{1}(visible),'--', ...
        'Color',colors(condition,:),'LineWidth',2.2, ...
        'DisplayName','DACF: fixed \omega');
    plot(ax,x(visible),curves{2}(visible),'-', ...
        'Color',colors(condition,:),'LineWidth',3, ...
        'DisplayName','DACF: fitted \omega');
    plot(ax,x(visible),curves{3}(visible),':', ...
        'Color',[0.12,0.12,0.12],'LineWidth',2.8, ...
        'DisplayName','MSD fit (method c)');
    xlim(ax,[0,40]);
    ylim(ax,[0,1.12*max([curves{1}(visible), ...
        curves{2}(visible),curves{3}(visible)])]);
    title(ax,condition_names(condition),'FontWeight','normal');
    xlabel(ax,'Modeled D_{eff} (\mum^2/s)','Interpreter','tex');
    if condition == 1
        ylabel(ax,'Probability density (s/\mum^2)','Interpreter','tex');
        legend(ax,'Location','northeast','Box','off','FontSize',12);
    end
    set(ax,'FontName','Arial','FontSize',14,'LineWidth',1.5, ...
        'TickDir','out','Box','on');
    ax.Toolbar.Visible = 'off';
    fprintf('%s PDF modes: fixed DACF %.3f; joint DACF %.3f; MSD fit %.3f um^2/s\n', ...
        condition_names(condition),modes);
end
title(layout,'SG3/11 trajectories: jointly fitting D_r and \omega to the DACF', ...
    'FontName','Arial','FontSize',17,'FontWeight','normal', ...
    'Interpreter','tex');
drawnow;
output_stem = fullfile(base_dir,'pico_diffusion_sg3_11_joint_dacf');
exportgraphics(fig,output_stem + ".png",'Resolution',400);
exportgraphics(fig,output_stem + ".pdf",'ContentType','vector');
savefig(fig,output_stem + ".fig");

function [t, correlation, pairs, omega] = track_dacf(xy,dt)
    displacement = diff(xy,1,1);
    lengths = hypot(displacement(:,1),displacement(:,2));
    valid = isfinite(lengths) & lengths > 0.05;
    turns_valid = valid(1:end-1) & valid(2:end);
    first = displacement(1:end-1,:);
    second = displacement(2:end,:);
    turns = atan2(first(:,1).*second(:,2)-first(:,2).*second(:,1), ...
        sum(first.*second,2));
    omega = abs(mean(turns(turns_valid))/dt);
    direction = nan(size(displacement));
    direction(valid,:) = displacement(valid,:)./lengths(valid);
    n = size(direction,1);
    t = (1:60)'*dt;
    correlation = nan(60,1);
    pairs = zeros(60,1);
    for lag = 1:60
        a = direction(1:n-lag,:);
        b = direction(1+lag:n,:);
        valid_pair = all(isfinite(a),2) & all(isfinite(b),2);
        pairs(lag) = sum(valid_pair);
        if pairs(lag) > 0
            correlation(lag) = mean(sum( ...
                a(valid_pair,:).*b(valid_pair,:),2));
        end
    end
    valid_lags = isfinite(correlation) & pairs>0;
    t = t(valid_lags);
    correlation = correlation(valid_lags);
    pairs = pairs(valid_lags);
end

function [Dr, amplitude] = fit_dr(t,y,weights,omega,free_amplitude)
    search = unique([0,linspace(0,75,100), ...
        logspace(-6,log10(75),120)]);
    scores = arrayfun(@(d) score_dr(d,t,y,weights,omega, ...
        free_amplitude),search);
    [~,best] = min(scores);
    low = search(max(1,best-1));
    high = search(min(numel(search),best+1));
    options = optimset('Display','off','TolX',1e-10);
    candidate = fminbnd(@(d) score_dr(d,t,y,weights,omega, ...
        free_amplitude),low,high,options);
    choices = [0,search(best),candidate,75];
    [~,choice] = min(arrayfun(@(d) score_dr(d,t,y,weights, ...
        omega,free_amplitude),choices));
    Dr = choices(choice);
    [~,amplitude] = score_dr(Dr,t,y,weights,omega,free_amplitude);
end

function [score, amplitude] = score_dr( ...
    Dr,t,y,weights,omega,free_amplitude)
    base = exp(-Dr*t).*cos(omega*t);
    if free_amplitude
        amplitude = min(1.2,max(0,sum(weights.*y.*base)/ ...
            sum(weights.*base.^2)));
    else
        amplitude = 1;
    end
    score = sum(weights.*(y-amplitude*base).^2);
end

function [Dr,omega] = fit_dr_free_omega( ...
    t,y,weights,prior_Dr,prior_omega,ensemble_omega)
    omega_bound = pi/(1/15);
    seeds = [prior_Dr,prior_omega; prior_Dr,ensemble_omega; ...
        max(0.1,prior_Dr),prior_omega*0.8; ...
        max(0.1,prior_Dr),prior_omega*1.2; ...
        0.1,ensemble_omega];
    options = optimoptions('lsqnonlin','Display','off', ...
        'MaxFunctionEvaluations',800);
    best_score = Inf;
    best = [NaN,NaN];
    for j = 1:size(seeds,1)
        seed = min(max(seeds(j,:),[0,0]),[75,omega_bound]);
        objective = @(q) sqrt(weights).* ...
            (y-exp(-q(1)*t).*cos(q(2)*t));
        [q,~,residual] = lsqnonlin(objective,seed,[0,0], ...
            [75,omega_bound],options);
        score = sum(residual.^2);
        if score < best_score
            best_score = score;
            best = q;
        end
    end
    Dr = best(1);
    omega = best(2);
end

function r2 = fit_r2(t,y,Dr,omega,amplitude)
    predicted = amplitude*exp(-Dr*t).*cos(omega*t);
    r2 = 1-sum((y-predicted).^2)/sum((y-mean(y)).^2);
end
