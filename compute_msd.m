function [D, msd, tau] = compute_msd(x, y, dt)
% Compute MSD and extract diffusion coefficient
% Inputs: x, y position vectors, dt timestep

    n    = length(x);
    lags = 1 : floor(n * 0.2);   % use up to 20% of trajectory
    msd  = zeros(size(lags));

    for k = 1:length(lags)
        lag     = lags(k);
        dx      = x(1+lag:end) - x(1:end-lag);
        dy      = y(1+lag:end) - y(1:end-lag);
        msd(k)  = mean(dx.^2 + dy.^2);
    end

    tau = lags * dt;

    % Fit line to long-time regime (last 50% of computed lags)
    % MSD = 4*D*tau  ->  slope = 4D
    long_idx = floor(length(tau)*0.5) : length(tau);
    p        = polyfit(tau(long_idx), msd(long_idx), 1);
    D        = p(1) / 4;

    fprintf('D_eff = %.4f um^2/s\n', D);

    % Plot
    figure; hold on; box on;
    plot(tau, msd, 'Color', [0.2 0.4 0.8], 'LineWidth', 1.5);
    plot(tau(long_idx), polyval(p, tau(long_idx)), ...
        'r--', 'LineWidth', 1.5);

    % Slope guide lines for reference
    tau_ref = tau(round(end/2));
    msd_ref = msd(round(end/2));
    plot(tau, msd_ref*(tau/tau_ref).^1, 'k:', 'LineWidth', 1);  % α=1
    plot(tau, msd_ref*(tau/tau_ref).^2, 'k--', 'LineWidth', 1); % α=2

    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('\tau (s)'); ylabel('MSD (\mum^2)');
    legend({'MSD', sprintf('fit: D=%.3f \\mum^2/s', D), ...
        '\alpha=1 (diffusive)', '\alpha=2 (ballistic)'}, ...
        'Location', 'northwest');
    title('Mean squared displacement');
end