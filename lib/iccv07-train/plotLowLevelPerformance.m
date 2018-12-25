function plotLowLevelPerformance

load('./results/classifierComparisonResult_stage3.mat');

figure(1), hold off
plot(pr_pb.r, pr_pb.p, 'b:', 'LineWidth', 4);

hold on
plot(pr_im.r, pr_im.p, 'g--', 'LineWidth', 3);
plot(pr_all.r, pr_all.p, 'r', 'LineWidth', 3);

legend({'Pb Only', 'Pb+Edge/Region Cues', 'Pb+Edge/Region+3D Cues'}, ...
    'Location', 'NorthEast', 'FontSize', 12);
title('Boundary Classification', 'FontSize', 20);
xlabel('Recall', 'FontSize', 18);
ylabel('Precision', 'FontSize', 18);
axis([0 1 0 1])
set(gca, 'FontSize', 14);