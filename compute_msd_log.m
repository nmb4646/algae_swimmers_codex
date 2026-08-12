function [D, alpha, msd, tau] = compute_msd_log(x, y, dt, Tm, wfrac)
% Compute MSD and extract diffusion coefficient
% Inputs: x, y position vectors, dt timestep
% Outputs: D     — effective diffusion coefficient (um^2/s)
%          alpha — scaling exponent (should be ~1 for diffusion)
%          msd   — mean squared displacement
%          tau   — lag times (s)

    n    = length(x);
    lags = 1 : floor(wfrac*n);
    msd  = zeros(size(lags));

    for k = 1:length(lags)
        lag    = lags(k);
        dx     = x(1+lag:end) - x(1:end-lag);
        dy     = y(1+lag:end) - y(1:end-lag);
        msd(k) = mean(dx.^2 + dy.^2);
    end

    tau = lags * dt;

    % Fit in log-log space over last 50% of computed lags
    % log(MSD) = alpha*log(tau) + log(4D)

    %long_idx = floor(length(tau)*0.5) : length(tau);
    long_idx = floor(Tm*.1):floor(.7*Tm);
    p        = polyfit(log(tau(long_idx)), log(msd(long_idx)), 1);
    alpha    = p(1);
    D        = exp(p(2)) / 4;

    fprintf('Alpha = %.2f\n', alpha);
    fprintf('D_eff = %.4e um^2/s\n', D);
   

    if abs(alpha - 1) > 0.3
        warning('Alpha = %.2f is far from 1 — fit may be unreliable', alpha);
    end
    if false
        % Plot
        figure; hold on; box on;
        scatter(tau, msd, 'Color', [0.2 0.4 0.8], 'LineWidth', 1.5);
        plot(tau(long_idx), exp(polyval(p, log(tau(long_idx)))), ...
            'r--', 'LineWidth', 1.5);
    
        % Slope guide lines for reference
        tau_ref = tau(round(end/2));
        msd_ref = msd(round(end/2));
        plot(tau, msd_ref*(tau/tau_ref).^1, 'k:',  'LineWidth', 1);
        plot(tau, msd_ref*(tau/tau_ref).^2, 'k--', 'LineWidth', 1);
    
        set(gca, 'XScale', 'log', 'YScale', 'log');
        xlabel('\tau (s)'); ylabel('MSD (\mum^2)');
        legend({'MSD', sprintf('fit: \\alpha=%.2f, D=%.3f \\mum^2/s', alpha, D), ...
            '\alpha=1 (diffusive)', '\alpha=2 (ballistic)'}, ...
            'Location', 'northwest');
        title(sprintf('Mean squared displacement — \\alpha=%.2f, D=%.4f \\mum^2/s', alpha, D));
    end

end