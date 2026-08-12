clear; clc; close all;

%%% SETTINGS
dir_m = './data/algae_motion/mucus_micro/';
n = 13;                          % trajectory number to analyze

dt  = 1/16.5;                    % s
pix = 0.24e-6;                   % m / pixel

go = 3;                          % sgolay order
gw = 11;                         % sgolay window

% ---- trap detection parameters ----
win_sec      = .5;              % sliding window size [s]
min_trap_sec = .5;              % minimum trapped duration to keep [s]
max_gap_sec  = 0.75;             % fill small free gaps inside trapped segments [s]
do_smooth_xy = true;             % smooth trajectory before detection
% -----------------------------------

%%% LOAD TRAJECTORY
filename = dir_m + sprintf("Microalgae in mucus (%i).csv", n);

warning('off','all');
T = readtable(filename);
warning('on','all');

x = (T.X - T.X(1)) * pix;
y = (T.Y - T.Y(1)) * pix;

if numel(x) < gw
    error('Trajectory is shorter than Savitzky-Golay window.');
end

if do_smooth_xy
    x = sgolayfilt(x, go, gw);
    y = sgolayfilt(y, go, gw);
end

t = (0:numel(x)-1)' * dt;
N = numel(x);

%%% CHECK TRAJECTORY LENGTH
w = max(5, 2*floor((win_sec/dt)/2)+1);   % make odd window length
halfW = floor(w/2);

if N < w
    error('Trajectory is shorter than the trap-detection window.');
end

%%% SLIDING-WINDOW FEATURES
% trapped windows should show:
%   - small radius of gyration
%   - small net displacement
%   - low straightness (back-and-forth motion)

Rg           = nan(N,1);   % local radius of gyration
netDisp      = nan(N,1);   % end-to-end displacement across window
pathLen      = nan(N,1);   % total path length in window
straightness = nan(N,1);   % netDisp / pathLen

for i = 1+halfW : N-halfW
    idx = (i-halfW):(i+halfW);

    xx = x(idx);
    yy = y(idx);

    % local radius of gyration
    xc = mean(xx);
    yc = mean(yy);
    rr2 = (xx - xc).^2 + (yy - yc).^2;
    Rg(i) = sqrt(mean(rr2));

    % local path length
    dxw = diff(xx);
    dyw = diff(yy);
    stepLen = hypot(dxw, dyw);
    pathLen(i) = sum(stepLen);

    % local net displacement
    netDisp(i) = hypot(xx(end) - xx(1), yy(end) - yy(1));

    % straightness
    straightness(i) = netDisp(i) / max(pathLen(i), eps);
end

valid = isfinite(Rg) & isfinite(netDisp) & isfinite(pathLen) & isfinite(straightness);

%%% UNSUPERVISED 2-STATE CLASSIFICATION
% Features are log-transformed to compress dynamic range.
F = [log10(Rg(valid) + eps), ...
     log10(netDisp(valid) + eps), ...
     log10(straightness(valid) + eps)];

% z-score standardization
muF = mean(F,1);
sdF = std(F,[],1);
Z = (F - muF) ./ max(sdF, eps);

% k-means into 2 states
rng(1);
idxK = kmeans(Z, 2, 'Replicates', 20, 'Display', 'off');

% choose the more confined cluster as "trapped"
c1 = mean(F(idxK==1,:), 1);
c2 = mean(F(idxK==2,:), 1);

score1 = sum(c1);   % smaller => more confined
score2 = sum(c2);

if score1 < score2
    trapCluster = 1;
else
    trapCluster = 2;
end

isTrapped = false(N,1);
isTrapped(valid) = (idxK == trapCluster);

%%% CLEAN BINARY STATE
min_trap_frames = max(1, round(min_trap_sec / dt));
max_gap_frames  = max(0, round(max_gap_sec / dt));

% fill short free gaps inside trapped segments
isTrapped = fillShortFalseGaps(isTrapped, max_gap_frames);

% remove very short trapped blips
isTrapped = removeShortTrueRuns(isTrapped, min_trap_frames);

%%% FIND TRAPPED SEGMENTS
[segStart, segEnd] = logicalSegments(isTrapped);
segDur = (segEnd - segStart + 1) * dt;

% longest trapped segment
trapMain = false(N,1);
if ~isempty(segStart)
    [~, imax] = max(segDur);
    trapMain(segStart(imax):segEnd(imax)) = true;
end

%%% EXTRACT LONGEST TRAPPED SUBTRAJECTORY
if any(trapMain)
    x_trap = x(trapMain);
    y_trap = y(trapMain);
    t_trap = t(trapMain);

    fprintf('Longest trapped segment: frames %d to %d (%.2f s)\n', ...
        find(trapMain,1,'first'), find(trapMain,1,'last'), numel(t_trap)*dt);
else
    x_trap = [];
    y_trap = [];
    t_trap = [];
    fprintf('No trapped segment detected in this trajectory.\n');
end

%%% VISUALIZATION
figure('Color','w','Position',[100 80 1250 850]);
tl = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

% ---------------------------------------------------
% PANEL 1: trajectory in xy
% ---------------------------------------------------
nexttile(tl,[2 1]); hold on; box on; axis equal;

% plot full trajectory once
plot(x, y, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1.2);

% overlay each contiguous trapped segment in red
for k = 1:numel(segStart)
    idx = segStart(k):segEnd(k);
    plot(x(idx), y(idx), 'r-', 'LineWidth', 2.0);
end

% highlight longest trapped segment a bit thicker
if any(trapMain)
    idxMain = find(trapMain);
    plot(x(idxMain), y(idxMain), '-', 'Color', [0.85 0 0], 'LineWidth', 3.0);
end

plot(x(1), y(1), 'go', 'MarkerFaceColor','g', 'MarkerSize',7);
plot(x(end), y(end), 'ko', 'MarkerFaceColor','k', 'MarkerSize',7);

xlabel('x [m]');
ylabel('y [m]');
title(sprintf('Trajectory %d: detected trapped regime', n));

legend_entries = {'full trajectory'};
legend_handles = plot(nan, nan, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1.2);

if ~isempty(segStart)
    h_trap = plot(nan, nan, 'r-', 'LineWidth', 2.0);
    legend_handles(end+1) = h_trap;
    legend_entries{end+1} = 'trapped segments';
end

h_start = plot(nan, nan, 'go', 'MarkerFaceColor','g', 'MarkerSize',7);
legend_handles(end+1) = h_start;
legend_entries{end+1} = 'start';

h_end = plot(nan, nan, 'ko', 'MarkerFaceColor','k', 'MarkerSize',7);
legend_handles(end+1) = h_end;
legend_entries{end+1} = 'end';

if any(trapMain)
    h_main = plot(nan, nan, '-', 'Color', [0.85 0 0], 'LineWidth', 3.0);
    legend_handles(end+1) = h_main;
    legend_entries{end+1} = 'longest trapped segment';
end

legend(legend_handles, legend_entries, 'Location', 'best');

% ---------------------------------------------------
% PANEL 2: local confinement features
% ---------------------------------------------------
nexttile; hold on; box on;

yyaxis left
yl1 = safeLimits(Rg);
ylim(yl1);
shadeSegments(segStart, segEnd, t, yl1);
h1 = plot(t, Rg, 'b-', 'LineWidth', 1.3);
ylabel('R_g [m]');

yyaxis right
yl2 = safeLimits(netDisp);
ylim(yl2);
h2 = plot(t, netDisp, '-', 'Color', [0.2 0.6 0.2], 'LineWidth', 1.3);
ylabel('net displacement [m]');

xlabel('time [s]');
title('Sliding-window confinement features');

% ---------------------------------------------------
% PANEL 3: straightness
% ---------------------------------------------------
nexttile; hold on; box on;

yl3 = safeLimits(straightness);
ylim(yl3);
shadeSegments(segStart, segEnd, t, yl3);
h3 = plot(t, straightness, 'k-', 'LineWidth', 1.3);

xlabel('time [s]');
ylabel('straightness = net/path');
title('Low straightness is consistent with trapping');

% ---------------------------------------------------
% PANEL 4: binary state vs time
% ---------------------------------------------------
nexttile([1 2]); hold on; box on;

stairs(t, double(isTrapped), 'r-', 'LineWidth', 2);
ylim([-0.1 1.1]);
yticks([0 1]);
yticklabels({'free','trapped'});
xlabel('time [s]');
ylabel('state');
title('Detected trapped state');

for k = 1:numel(segStart)
    text(t(segStart(k)), 1.03, sprintf('%.2f s', segDur(k)), ...
        'Color', [0.7 0 0], 'FontSize', 9, ...
        'VerticalAlignment', 'bottom');
end

sgtitle(sprintf('Trap detection for Microalgae in mucus (%d)', n), ...
    'FontWeight', 'bold');

%%% OPTIONAL: print summary
fprintf('\nDetected trapped segments:\n');
if isempty(segStart)
    fprintf('  none\n');
else
    for k = 1:numel(segStart)
        fprintf('  segment %d: frames %d-%d, duration = %.2f s\n', ...
            k, segStart(k), segEnd(k), segDur(k));
    end
end

%%% =========================
% LOCAL FUNCTIONS
% =========================

function [s,e] = logicalSegments(state)
    state = logical(state(:));
    d = diff([false; state; false]);
    s = find(d == 1);
    e = find(d == -1) - 1;
end

function state = removeShortTrueRuns(state, minLen)
    if isempty(state) || minLen <= 1
        return;
    end

    [s,e] = logicalSegments(state);
    for ii = 1:numel(s)
        if (e(ii) - s(ii) + 1) < minLen
            state(s(ii):e(ii)) = false;
        end
    end
end

function state = fillShortFalseGaps(state, maxGap)
    if isempty(state) || maxGap <= 0
        return;
    end

    state = logical(state(:));

    % find runs of false values
    falseState = ~state;
    [s,e] = logicalSegments(falseState);

    for ii = 1:numel(s)
        gapLen = e(ii) - s(ii) + 1;

        leftIsTrap  = (s(ii) > 1) && state(s(ii)-1);
        rightIsTrap = (e(ii) < numel(state)) && state(e(ii)+1);

        if gapLen <= maxGap && leftIsTrap && rightIsTrap
            state(s(ii):e(ii)) = true;
        end
    end
end

function yl = safeLimits(v)
    v = v(isfinite(v));
    if isempty(v)
        yl = [0 1];
        return;
    end

    vmin = min(v);
    vmax = max(v);

    if vmin == vmax
        pad = max(abs(vmin)*0.05, eps);
        yl = [vmin-pad, vmax+pad];
    else
        pad = 0.05 * (vmax - vmin);
        yl = [vmin-pad, vmax+pad];
    end
end

function shadeSegments(segStart, segEnd, t, yl)
    if isempty(segStart)
        return;
    end

    for kk = 1:numel(segStart)
        xs = [t(segStart(kk)) t(segEnd(kk)) t(segEnd(kk)) t(segStart(kk))];
        ys = [yl(1) yl(1) yl(2) yl(2)];
        patch(xs, ys, [1 0.85 0.85], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.35);
    end
end