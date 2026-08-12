% File name
close all; clear; clc;

blue = [.06,.43,.74];

dir_w = './data/algae_motion/buffer_micro/';
dir_m = './data/algae_motion/mucus_micro/';
dt = 1/16.5; pix = .24e-6;
Nw = 20;Nm=20;
% Scalar features to extract

mean_speed_w = zeros(Nw,1);
mean_curvature_w = zeros(Nw,1);
winding_number_w = zeros(Nw,1);
winding_rate_w = zeros(Nw,1);

loop_circularity_w = [];
loop_length_w = [];
loop_chirality_w = [];

mean_speed_m= zeros(Nm,1);
mean_curvature_m = zeros(Nm,1);
winding_number_m = zeros(Nm,1);
winding_rate_m = zeros(Nm,1);

loop_circularity_m = [];
loop_length_m = [];
loop_chirality_m = [];


% Vector features to extract
mean_velocity_w = zeros(Nw,2);
mean_velocity_m = zeros(Nm,2);
loop_velocity_w = [];
loop_velocity_m = [];

%Function features to extract

vacf_w = cell(Nw,1);
vacf_m = cell(Nm,1);



%Water feature extraction
for n = 1:Nm
    filename = dir_w + sprintf("Microalgae in water (%i).csv",n);
    warning('off','all');
    T = readtable(filename);
    warning('on','all');
    % Extract X and Y columns (make sure names match exactly)
    x = T.X; x = x-x(1);x=x*pix;
    y = T.Y; y = y-y(1);y=y*pix;
    % Combine into a 2xN matrix
    A = [x, y];   % transpose to make row vectors
    disp(n);disp(length(x))
    %smoothing
    go = 3;gw=11;
    mm = 1; x_unsmoothed=x;y_unsmoothed=y;
    x = sgolayfilt(x,go,gw); y = sgolayfilt(y,go,gw);x=movmean(x,mm);y=movmean(y,mm);
    dx = gradient(x)/dt;dy = gradient(y)/dt;
    ddx = gradient(dx)/dt;ddy = gradient(dy)/dt;
    speed = sqrt(dot([dx,dy],[dx,dy],2));
    
    k = movmean((dx .* ddy - dy .* ddx) ./ (dx.^2 + dy.^2).^(3/2),1);
    k(isnan(k))=0;
    theta = atan2(dy, dx);
    % unwrap to avoid jumps at ±pi
    theta_unwrapped = unwrap(theta);
    
    % total turning angle
    total_turn = theta_unwrapped(end) - theta_unwrapped(1);
    
    % winding number
    wn = total_turn / (2*pi);

    % report trajectory-level values
    mean_speed_w(n) = mean(speed);
    mean_curvature_w(n) = 1/mean(k);
    winding_number_w(n) = wn;
    winding_rate_w(n) = wn/length(x);
    mean_velocity_w(n,:) = mean([dx,dy],1);

    % loop segmentation

    theta_rel = theta_unwrapped - theta_unwrapped(1);
    total_rot = abs(theta_rel(end)) / (2*pi);
    n_full    = floor(total_rot);
    
    % --- Find first crossing of each 2π multiple ---
    crossings    = zeros(1, n_full + 1);
    crossings(1) = 1;
    for k = 1:n_full
        crossings(k+1) = find(abs(theta_rel) >= k*2*pi, 1, 'first');
    end
    % n_full complete loops, ignoring everything after the last crossing
    % boundaries has n_full+1 entries defining n_full segments
    boundaries = crossings;
    n_segments = n_full;   % strictly complete loops only
    
    for i = 1:n_segments
        idx = boundaries(i) : boundaries(i+1);
        xi  = x(idx);
        yi  = y(idx);
    
        % Arc coverage sanity check
        arc = abs(theta_rel(idx(end)) - theta_rel(idx(1))) / (2*pi);
        loop_chirality_w = [loop_chirality_w,sign(theta_rel(idx(end)) - theta_rel(idx(1)))];
        loop_length_w = [loop_length_w;sum(sqrt(dx(idx).^2 + dy(idx).^2))];
    
        % Circle fit
        [circul, r___] = fit_circle(xi, yi);
        loop_circularity_w = [loop_circularity_w;circul];

        % loop drift
        loop_velocity_w = [loop_velocity_w; [x(idx(end)) - x(idx(1)), y(idx(end)) - y(idx(1))]];
    end

    vacf_w{n} = vacf(dx,dy,dt);

end

%Mucus feature extraction
for n = 1:[]
        filename = dir_m + sprintf("Microalgae in mucus (%i).csv",n);
    warning('off','all');
    T = readtable(filename);
    warning('on','all');
    % Extract X and Y columns (make sure names match exactly)
    x = T.X; x = x-x(1);x=x*pix;
    y = T.Y; y = y-y(1);y=y*pix;

    %plot(x,y);figure;
    % Combine into a 2xN matrix
    A = [x, y];   % transpose to make row vectors
    %smoothing
    go = 3;gw=11;
    mm = 1;
    x = sgolayfilt(x,go,gw); y = sgolayfilt(y,go,gw);x=movmean(x,mm);y=movmean(y,mm);
    dx = gradient(x)/dt;dy = gradient(y)/dt;
    ddx = gradient(dx)/dt;ddy = gradient(dy)/dt;
    speed = sqrt(dot([dx,dy],[dx,dy],2));
    
    k = movmean((dx .* ddy - dy .* ddx) ./ (dx.^2 + dy.^2).^(3/2),1);
    k(isnan(k))=0;
    %disp("Mean radius: " + sprintf("%.2f",1/mean(k)))
    theta = atan2(dy, dx);
    % unwrap to avoid jumps at ±pi
    theta_unwrapped = unwrap(theta);
    
    % total turning angle
    total_turn = theta_unwrapped(end) - theta_unwrapped(1);
    
    % winding number
    wn = total_turn / (2*pi);

    %Report trajectory-level values
    mean_speed_m(n) = mean(speed);
    mean_curvature_m(n) = 1/mean(k);
    winding_number_m(n) = wn;
    winding_rate_m(n) = wn/length(x);
    mean_velocity_m(n,:) = mean([dx,dy],1);

    
    % loop segmentation

    theta_rel = theta_unwrapped - theta_unwrapped(1);
    total_rot = abs(theta_rel(end)) / (2*pi);
    n_full    = floor(total_rot);
    
    % --- Find first crossing of each 2π multiple ---
    crossings    = zeros(1, n_full + 1);
    crossings(1) = 1;
    for k = 1:n_full
        crossings(k+1) = find(abs(theta_rel) >= k*2*pi, 1, 'first');
    end
    % n_full complete loops, ignoring everything after the last crossing
    % boundaries has n_full+1 entries defining n_full segments
    boundaries = crossings;
    n_segments = n_full;   % strictly complete loops only
    
    for i = 1:n_segments
        idx = boundaries(i) : boundaries(i+1);
        xi  = x(idx);
        yi  = y(idx);
    
        % Arc coverage sanity check
        arc = abs(theta_rel(idx(end)) - theta_rel(idx(1))) / (2*pi);
        loop_chirality_m = [loop_chirality_m,sign(theta_rel(idx(end)) - theta_rel(idx(1)))];
        loop_length_m = [loop_length_m;sum(sqrt(dx(idx).^2 + dy(idx).^2))];
    
        % Circle fit
        [circul, r___] = fit_circle(xi, yi);
        loop_circularity_m = [loop_circularity_m;circul];

        % loop drift
        loop_velocity_m = [loop_velocity_m; [x(idx(end)) - x(idx(1)), y(idx(end)) - y(idx(1))]];
    end
    
    vacf_m{n} = vacf(dx,dy,dt);

end


zw = 0*winding_rate_w;zm = 0*winding_rate_m; 

if false % filter loop-level features to nearly circular loops
    th = .91;
    loop_circularity_w = loop_circularity_w(loop_circularity_w>th);
    loop_length_w = loop_length_w(loop_circularity_w>th);
    loop_chirality_w =loop_chirality_w(loop_circularity_w>th);
    loop_velocity_w = loop_velocity_w(loop_circularity_w>th,:);

    loop_circularity_m = loop_circularity_m(loop_circularity_m>th);
    loop_length_m = loop_length_m(loop_circularity_m>th);
    loop_chirality_m =loop_chirality_m(loop_circularity_m>th);
    loop_velocity_m = loop_velocity_m(loop_circularity_m>th,:);
end
zwl = zeros(length(loop_length_w),1); zml = zeros(length(loop_length_m),1);

%Rayleigh test for orientations/velocities

for orientation_calc = []
figure;tiledlayout(1,2);nexttile; plot_rose(loop_velocity_w);title(sprintf('Water, p = %.3f',prw)); nexttile;plot_rose(loop_velocity_m);title(sprintf('Mucus, p = %.3f',prm))
set(gcf,"Position",[100 100 900 600]);exportgraphics(gcf,'meandrift.png')
end


%set(gca,"XScale","linear")
x=x-x(1);y=y-y(1);



%compute_msd_log(xd,yd,dt)

figure;hold on;
%plot(x(boundaries),y(boundaries))


plot(x,y,LineWidth=2);scatter(x(1),y(1)); hold on;

title("Duration = " + dt*length(x))







for vacf_plot =[]
for n = []
    if n == 36
    plot(vacf_w{n},LineWidth=3,Color=[.2,.4,.8])
    else
    plot(vacf_w{n},LineWidth=2,Color=[.5,.5,.5,.2])
    end
    xlim([0,800])
    title('VACFs in water');xlabel('t');ylabel('Correlation')
end

set(gcf,"Position",[100 100 1000 800])
exportgraphics(gcf,'VACFs.png')
end

for mannwhitney = 1
[p,h,stats] = ranksum(mean_speed_w,mean_speed_m);
n1 = length(zwl);
n2 = length(zml);
U = stats.ranksum - n1*(n1+1)/2;
r = (2*U)/(n1*n2) - 1;

fprintf('ranksum = %.1f\n', stats.ranksum)
fprintf('U       = %.1f\n', U)
fprintf('r       = %.3f\n', r)
fprintf('p       = %.4f\n', p)
end

for smoothing =[]
%plot(x_unsmoothed,y_unsmoothed,LineWidth=2,color=[0.8500 0.3250 0.0980]);plot(x,y,LineWidth=2,Color=blue); hold on; scatter(x(1),y(1),70,blue,LineWidth=2,MarkerFaceColor='w')
plot(x_unsmoothed, y_unsmoothed, '-', ...
    'Color', [0.7 0.7 0.7 ,1], ...
    'LineWidth', 3);
% plot(x_unsmoothed, y_unsmoothed, '-', ...
%     'Color', [0.3 0.3 0.3 0.3], ...
%     'LineWidth', 4);
hold on;
% Smooth — bold, colored
plot(x, y, '-', ...
    'Color', [0.9,0.1,0.1], ...
    'LineWidth', 2);
axis off;
axis equal;
exportgraphics(gcf,'smoothing.png')
end

for osculating=[]
if false
    plot(x,y,LineWidth=1.5);hold on; scatter(x(1),y(1),30,[0, 0.4470, 0.7410]);axis equal;%figure;
    %velocities 
    %figure;
    %plot(k);ylim([-.1,.1])
end


if false
    %PLOT OSCULATING CIRCLE
    idx = 51; kappa=k;
    % Unit tangent
    T = [dx(idx); dy(idx)];
    T = T / norm(T);
    % Unit normal (90° CCW)
    N = [-T(2); T(1)];
    % Radius
    R = 1 / abs(kappa(idx));
    % Center of osculating circle (includes sign of curvature)
    C = [x(idx); y(idx)] + (1 / kappa(idx)) * N;
    % Circle points
    theta = linspace(0, 2*pi, 200);
    xc = C(1) + R * cos(theta);
    yc = C(2) + R * sin(theta);
    % Plot on existing figure
    plot(xc, yc, 'r--', 'LineWidth', 1.5);
    plot(x(idx), y(idx), 'ko', 'MarkerFaceColor','k');
    plot(C(1), C(2), 'ro', 'MarkerFaceColor','r');
end
end

for loopsegviz =[]
    % --- Colormap ---
    cmap = hsv(n_segments);
    circularity = zeros(1, n_segments);
    radii       = zeros(1, n_segments);
    for i = 1:n_segments
        idx = boundaries(i) : boundaries(i+1);
        xi  = x(idx);
        yi  = y(idx);
    
        % Arc coverage sanity check
        arc = abs(theta_rel(idx(end)) - theta_rel(idx(1))) / (2*pi);
    
        % Circle fit
        [circularity(i), radii(i)] = fit_circle(xi, yi);
    
        % Plot segment
        plot(xi, yi, 'Color', cmap(i,:), 'LineWidth', 2.5);
    
        % Annotate with circularity at loop midpoint
        % mid = round(length(idx)/2);
        % text(xi(mid), yi(mid), sprintf('C=%.2f\nR=%.1f', circularity(i), radii(i)), ...
        %     'FontSize', 7, 'Color', cmap(i,:)*0.7, 'HorizontalAlignment', 'center');
    end
    
    % Start marker only — no end marker since we're ignoring the tail
    plot(x(1), y(1), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
    
    colormap(cmap);
    cb = colorbar(fontsize=20);
    clim([0, n_segments-1]);
    cb.Label.String = 'Loop number';
    title(sprintf('%d full loops of %.2f total rotations (%.2f ignored)', ...
        n_full, total_rot, total_rot - n_full),FontSize=20);
    xlabel('x',FontSize=20); ylabel('y','FontSize',20);
end
set(gcf,'Position',[1000,1000,1000,800])
axis equal
exportgraphics(gcf,'segviz1.png')

function [circularity, R, center] = fit_circle(xi, yi)
    A  = [xi(:), yi(:), ones(length(xi), 1)];
    b  = xi(:).^2 + yi(:).^2;
    c  = A \ b;
    xc = c(1)/2;
    yc = c(2)/2;
    R  = sqrt(c(3) + xc^2 + yc^2);
    center = [xc, yc];

    dist      = sqrt((xi(:)-xc).^2 + (yi(:)-yc).^2);
    residuals = dist - R;
    circularity = 1 - (std(residuals) / R);
end

function [vacf_out,tau_d_out] = vacf(dx,dy,dt)
    speed = sqrt(dx.^2 + dy.^2);
    vx_hat = dx ./ speed;
    vy_hat = dy ./ speed;
    
    % Autocorrelation of each component, average them
    [cx, lags] = xcorr(vx_hat, 'normalized');
    [cy, ~]    = xcorr(vy_hat, 'normalized');
    
    pos  = lags >= 0;
    C    = (cx(pos) + cy(pos)) / 2;
    tau  = lags(pos) * dt;
    
    % --- Extract persistence time: first zero crossing ---
    sign_changes = find(diff(sign(C)));
    if ~isempty(sign_changes)
        % Linear interpolation between the two points straddling zero
        i1 = sign_changes(1);
        i2 = i1 + 1;
        tau_p = tau(i1) + (tau(i2)-tau(i1)) * (-C(i1))/(C(i2)-C(i1));
        %fprintf('Persistence time tau_p = %.4f s\n', tau_p);
    else
        tau_p = NaN;
        %fprintf('No zero crossing found — trajectory may be too short\n');
    end
    
    % --- Extract loop period: first peak after zero crossing ---
    if ~isnan(tau_p)
        after_zero = find(tau > tau_p);
        if ~isempty(after_zero)
            C_after   = C(after_zero);
            tau_after = tau(after_zero);
            [~, peak_idx] = max(C_after);
            loop_period = tau_after(peak_idx);
            %fprintf('Loop period = %.4f s\n', loop_period);
        else
            loop_period = NaN;
        end
    end
    
    
    % --- Plot ---
    if false
        figure; hold on; box on;
        plot(tau, C, 'Color', [0.2 0.4 0.8], 'LineWidth', 1.5);
        yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
        if ~isnan(tau_p)
            xline(tau_p, 'r--', 'LineWidth', 1, 'Label', '\tau_p');
        end
        xlabel('Lag \tau (s)');
        ylabel('C(\tau)');
        title('Velocity autocorrelation function');
        xlim([0, tau(end)]);
        ylim([-1.1, 1.1]);
    end
    
    vacf_out = C;

end

function [p, R, mean_angle] = rayleigh_test(velocities)
% RAYLEIGH_TEST  Tests whether 2D velocity vectors have a preferred direction
%
% Inputs:
%   velocities  — Nx2 array of [vx, vy] velocities
%
% Outputs:
%   p           — p-value (small = significant preferred direction)
%   R           — mean resultant length (0 = uniform, 1 = all identical)
%   mean_angle  — mean direction in radians

angles = atan2(velocities(:,2), velocities(:,1));

n     = length(angles);
R     = abs(sum(exp(1i * angles))) / n;
z     = n * R^2;
p     = exp(-z) * (1 + (2*z - z^2)/(4*n) - (24*z - 132*z^2 + 76*z^3 - 9*z^4)/(288*n^2));

mean_angle = angle(sum(exp(1i * angles)));

fprintf('Mean resultant length R = %.3f\n', R);
fprintf('Mean drift direction    = %.1f deg\n', rad2deg(mean_angle));
fprintf('Rayleigh test p         = %.4f\n', p);

end

function plot_rose(velocities, titlestr)
% PLOT_ROSE  Rose plot of 2D velocity directions with mean resultant vector
%
% Inputs:
%   velocities  — Nx2 array of [vx, vy] velocities
%   titlestr    — optional title string

if nargin < 2
    titlestr = 'Drift directions';
end



angles     = atan2(velocities(:,2), velocities(:,1));
n          = length(angles);
R          = abs(sum(exp(1i * angles))) / n;
mean_angle = angle(sum(exp(1i * angles)));


polarhistogram(angles, 12, ...
    'FaceColor',      [0.2 0.4 0.8], ...
    'FaceAlpha',       0.6, ...
    'EdgeColor',      [0.2 0.4 0.8] * 0.7);

hold on;

% Mean resultant vector — length scaled to plot radius
ax    = gca;
r_max = max(ax.RLim);
polarplot([mean_angle, mean_angle], [0, R * r_max], ...
    'r-', 'LineWidth', 2.5);
polarplot(mean_angle, R * r_max, ...
    'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);

end