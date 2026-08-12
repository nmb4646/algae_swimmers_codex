close all; clear; clc;

tau = logspace(-3,3,1000);
f = (tau.^2)./(tau+1);

figure; hold on
loglog(tau, f, 'LineWidth', 3, 'HandleVisibility','off')

% dashed asymptotic fit lines
loglog(tau, tau.^2, '--', 'LineWidth', 1.5, 'DisplayName', 'Ballistic transport, \alpha = 2')
loglog(tau, tau, '--', 'LineWidth', 1.5, 'DisplayName', 'Diffusive transport, \alpha = 1')
xlabel('\tau')
ylabel('<r^2>')
legend('Location','best')
title("Example MSD")

set(gca,"XScale","log")
set(gca,"YScale","log")
set(gca,"FontSize",15)
