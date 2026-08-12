vals = [6.05e-12, 6.3e-11, 1.05e-10];

figure;
b=bar(vals);
set(gca, 'YScale', 'log');

xticks(1:3);
xticklabels({'0%','1%','2%'});
xlabel("Mucin percentage")
ylabel('D_{eff} (m^2/s)');
title('Analytical Diffusion Coefficients');

ylim([9e-13,1e-9])

b.FaceColor = 'flat';
b.CData = [
    0.6 0.4 0.8
    0.9 0.3 0.2
    0.3 0.7 0.3
];
set(gca,"FontSize",18)

grid on;