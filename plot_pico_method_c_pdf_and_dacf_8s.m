clear; clc; close all;

% Standalone examples: SG3/11 method-c D_eff PDF with swimmer-bootstrap
% bands, and the original ensemble DACF recomputed to twice its 4-s lag.
base_dir = fileparts(mfilename('fullpath'));
diffusion_source = fullfile(base_dir, ...
    'compare_pico_diffusion_methods_abc_sg1_3_results.mat');
dacf_source = fullfile(base_dir, 'analysis_2jklm_msd_fit_results.mat');
condition_dirs = ["picoalgae 0% mucin"; "picoalgae 2% mucin"];
condition_names = ["0% mucin"; "2% mucin"];
colors = [0.12,0.47,0.71; 1.00,0.25,0.39];
dt = 1/15;
pixel_size_um = 0.48;
bootstrap_samples = 1000;
bootstrap_seed = 20260928;
pdf_xlimits = [0,40]; % Peak-region view, without renormalization.
dacf_max_lag_s = 8;
min_direction_step_um = 0.05;

%%%% METHOD-C PDF: SAME SG3/11 COHORT, KDE GRID, AND BANDWIDTH
S = load(diffusion_source, 'settings','group_index','c_values','x_grid');
variant = find(S.settings.sg_order==3 & ...
    S.settings.sg_window_frames==11);
assert(isscalar(variant), 'SG3/11 is missing from the saved comparison.');
grid = S.x_grid(:)';
visible = grid>=pdf_xlimits(1) & grid<=pdf_xlimits(2);
x = grid(visible);
h = S.settings.kde_bandwidth_log;
pdf_mean = cell(2,1);
pdf_sem = cell(2,1);
pdf_lower = cell(2,1);
pdf_upper = cell(2,1);

old_rng = rng;
restore_rng = onCleanup(@() rng(old_rng));
for condition = 1:2
    values = S.c_values(variant,S.group_index{condition});
    values = values(:);
    assert(all(isfinite(values) & values>0));
    kernels = zeros(numel(values),numel(x));
    positive = x>0;
    z = (log(x(positive))-log(values))/h;
    kernels(:,positive) = exp(-0.5*z.^2) ./ ...
        (sqrt(2*pi)*h*x(positive));
    pdf_mean{condition} = mean(kernels,1);
    rng(bootstrap_seed+condition-1,'twister');
    n = numel(values);
    pdf_sem{condition} = std(kernels,0,1)/sqrt(n);
    draws = randi(n,n,bootstrap_samples);
    columns = repelem((1:bootstrap_samples)',n);
    weights = accumarray([draws(:),columns],1,[n,bootstrap_samples])/n;
    sampled_pdfs = weights'*kernels;
    limits = prctile(sampled_pdfs,[2.5,97.5],1);
    pdf_lower{condition} = limits(1,:);
    pdf_upper{condition} = limits(2,:);
    fprintf('%s method-c PDF: %d swimmers; SG3/11; 95%% pointwise bootstrap band\n', ...
        condition_names(condition),n);
end

pdf_peak = nan(2,1);
for condition = 1:2
    [~,index] = max(pdf_mean{condition});
    pdf_peak(condition) = x(index);
end

fig_pdf = figure('Color','w','Position',[90,90,1050,620]);
ax = axes(fig_pdf); hold(ax,'on');
for condition = 1:2
    shade_band(ax,x,pdf_lower{condition},pdf_upper{condition}, ...
        colors(condition,:),0.18);
    xline(ax,pdf_peak(condition),'--','Color',colors(condition,:), ...
        'LineWidth',1.5,'HandleVisibility','off');
    plot(ax,x,pdf_mean{condition},'Color',colors(condition,:), ...
        'LineWidth',2.3,'DisplayName',sprintf('%s (peak %.1f)', ...
        condition_names(condition),pdf_peak(condition)));
end
xlim(ax,pdf_xlimits);
ylim(ax,[0,1.08*max([pdf_upper{1},pdf_upper{2}])]);
xlabel(ax,'Modeled D_{eff} (\mum^2/s)','Interpreter','tex');
ylabel(ax,'Probability density (s/\mum^2)','Interpreter','tex');
title(ax,'Method c: single-swimmer MSD fits','FontWeight','normal');
legend(ax,'Location','northeast','Box','off');
text(ax,0.97,0.80,'95% pointwise swimmer-bootstrap bands', ...
    'Units','normalized','HorizontalAlignment','right', ...
    'FontSize',11,'Color',[0.35,0.35,0.35]);
style_axes(ax);
pdf_stem = fullfile(base_dir,'pico_method_c_sg3_11_pdf_bootstrap');
exportgraphics(fig_pdf,pdf_stem + ".png",'Resolution',400);
exportgraphics(fig_pdf,pdf_stem + ".pdf",'ContentType','vector');
savefig(fig_pdf,pdf_stem + ".fig");

% A second version uses the same pointwise mean +/- one standard error
% across swimmer kernels, analogous to the curvature PDF panel.
fig_sem = figure('Color','w','Position',[90,90,1050,620]);
ax = axes(fig_sem); hold(ax,'on');
sem_upper = cell(2,1);
for condition = 1:2
    sem_lower = max(0,pdf_mean{condition}-pdf_sem{condition});
    sem_upper{condition} = pdf_mean{condition}+pdf_sem{condition};
    shade_band(ax,x,sem_lower,sem_upper{condition}, ...
        colors(condition,:),0.18);
    xline(ax,pdf_peak(condition),'--','Color',colors(condition,:), ...
        'LineWidth',1.5,'HandleVisibility','off');
    plot(ax,x,pdf_mean{condition},'Color',colors(condition,:), ...
        'LineWidth',2.3,'DisplayName',sprintf('%s (peak %.1f)', ...
        condition_names(condition),pdf_peak(condition)));
end
xlim(ax,pdf_xlimits);
ylim(ax,[0,1.08*max([sem_upper{1},sem_upper{2}])]);
xlabel(ax,'Modeled D_{eff} (\mum^2/s)','Interpreter','tex');
ylabel(ax,'Probability density (s/\mum^2)','Interpreter','tex');
title(ax,'Method c: single-swimmer MSD fits','FontWeight','normal');
legend(ax,'Location','northeast','Box','off');
text(ax,0.97,0.80,'Mean \pm 1 SEM across swimmers', ...
    'Units','normalized','HorizontalAlignment','right', ...
    'Interpreter','tex','FontSize',11,'Color',[0.35,0.35,0.35]);
style_axes(ax);
sem_stem = fullfile(base_dir,'pico_method_c_sg3_11_pdf_sem');
exportgraphics(fig_sem,sem_stem + ".png",'Resolution',400);
exportgraphics(fig_sem,sem_stem + ".pdf",'ContentType','vector');
savefig(fig_sem,sem_stem + ".fig");

%%%% ORIGINAL PANEL-L DACF: SAME SG2/5 SMOOTHING, EXTENDED TO 8 S
Q = load(dacf_source,'results');
max_lag = round(dacf_max_lag_s/dt);
lag_s = (0:max_lag)'*dt;
dacf_mean = cell(2,1);
dacf_sem = cell(2,1);
dacf_count = cell(2,1);
for condition = 1:2
    files = Q.results(condition).files;
    correlations = nan(max_lag+1,numel(files));
    for swimmer = 1:numel(files)
        path = fullfile(base_dir,'data','Picoalgae in 0-2% mucin', ...
            condition_dirs(condition),files(swimmer));
        xy = read_track_as_original(path,dt,pixel_size_um);
        xy = [sgolayfilt(xy(:,1),2,5),sgolayfilt(xy(:,2),2,5)];
        correlations(:,swimmer) = track_dacf( ...
            xy,max_lag,min_direction_step_um);
    end
    dacf_mean{condition} = mean(correlations,2,'omitnan');
    dacf_count{condition} = sum(isfinite(correlations),2);
    dacf_sem{condition} = std(correlations,0,2,'omitnan') ./ ...
        sqrt(dacf_count{condition});
    dacf_sem{condition}(dacf_count{condition}<2) = NaN;
    saved_mean = Q.results(condition).dacf_mean;
    saved_sem = Q.results(condition).dacf_sem;
    assert(max(abs(dacf_mean{condition}(1:numel(saved_mean))-saved_mean))<1e-9);
    assert(max(abs(dacf_sem{condition}(1:numel(saved_sem))-saved_sem))<1e-9);
    fprintf('%s DACF: %d swimmers at zero lag; %d at 8 s\n', ...
        condition_names(condition),dacf_count{condition}(1), ...
        dacf_count{condition}(end));
end

fig_dacf = figure('Color','w','Position',[90,90,1050,620]);
ax = axes(fig_dacf); hold(ax,'on');
all_bounds = 0;
for condition = 1:2
    lower = dacf_mean{condition}-dacf_sem{condition};
    upper = dacf_mean{condition}+dacf_sem{condition};
    shade_band(ax,lag_s,lower,upper,colors(condition,:),0.18);
    plot(ax,lag_s,dacf_mean{condition}, ...
        'Color',colors(condition,:),'LineWidth',2.3, ...
        'DisplayName',condition_names(condition));
    all_bounds = [all_bounds;lower(:);upper(:)]; %#ok<AGROW>
end
yline(ax,0,':','Color',[0.35,0.35,0.35], ...
    'LineWidth',0.9,'HandleVisibility','off');
xlim(ax,[0,dacf_max_lag_s]);
all_bounds = all_bounds(isfinite(all_bounds));
padding = 0.05*range(all_bounds);
ylim(ax,[min(-0.5,floor(10*(min(all_bounds)-padding))/10), ...
    max(1,ceil(10*(max(all_bounds)+padding))/10)]);
xlabel(ax,'\tau (s)','Interpreter','tex');
ylabel(ax,'Direction correlation');
title(ax,'Ensemble DACF through 8 s','FontWeight','normal');
legend(ax,'Location','northeast','Box','off');
style_axes(ax);
dacf_stem = fullfile(base_dir,'pico_ensemble_dacf_8s');
exportgraphics(fig_dacf,dacf_stem + ".png",'Resolution',400);
exportgraphics(fig_dacf,dacf_stem + ".pdf",'ContentType','vector');
savefig(fig_dacf,dacf_stem + ".fig");

function xy = read_track_as_original(path,dt,pixel_size_um)
    T = readtable(path,'VariableNamingRule','preserve');
    variables = string(T.Properties.VariableNames);
    if ismember("FrameNumber",variables)
        frames = double(T.FrameNumber(:));
    else
        frames = (0:height(T)-1)';
    end
    x = double(T.X(:));
    y = double(T.Y(:));
    valid = isfinite(frames) & isfinite(x) & isfinite(y);
    frames = frames(valid); x = x(valid); y = y(valid);
    [frames,order] = sort(frames);
    x = x(order); y = y(order);
    [frames,index] = unique(frames,'stable');
    x = x(index); y = y(index);
    t = (frames-frames(1))*dt;
    x = (x-x(1))*pixel_size_um;
    y = (y-y(1))*pixel_size_um;
    if any(abs(diff(t)-dt)>100*eps(max(t(end),dt)))
        uniform_t = (0:floor(t(end)/dt))'*dt;
        x = interp1(t,x,uniform_t,'pchip');
        y = interp1(t,y,uniform_t,'pchip');
    end
    xy = [x(:),y(:)];
end

function correlation = track_dacf(xy,max_lag,min_step_um)
    displacement = diff(xy,1,1);
    lengths = hypot(displacement(:,1),displacement(:,2));
    valid = isfinite(lengths) & lengths>min_step_um;
    direction = nan(size(displacement));
    direction(valid,:) = displacement(valid,:)./lengths(valid);
    n = size(direction,1);
    correlation = nan(max_lag+1,1);
    for lag = 0:min(max_lag,n-1)
        first = direction(1:n-lag,:);
        second = direction(1+lag:n,:);
        pair = all(isfinite(first),2) & all(isfinite(second),2);
        if any(pair)
            correlation(lag+1) = mean(sum( ...
                first(pair,:).*second(pair,:),2));
        end
    end
end

function shade_band(ax,x,lower,upper,color,alpha)
    good = isfinite(x) & isfinite(lower) & isfinite(upper);
    x = x(good); lower = lower(good); upper = upper(good);
    fill(ax,[x(:);flipud(x(:))], ...
        [lower(:);flipud(upper(:))],color, ...
        'FaceAlpha',alpha,'EdgeColor','none','HandleVisibility','off');
end

function style_axes(ax)
    ax.Position = [0.11,0.13,0.84,0.73];
    set(ax,'FontName','Arial','FontSize',14,'LineWidth',1.5, ...
        'TickDir','out','Box','off');
    ax.Toolbar.Visible = 'off';
end
