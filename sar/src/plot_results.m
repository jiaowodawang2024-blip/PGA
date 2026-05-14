function plot_results(raw, range_data, image, targets, params)
%PLOT_RESULTS Display raw, range-compressed, and focused SAR results.

figure('Name', 'Near-field MIMO SAR simulation', 'Color', 'w');

subplot(2, 2, 1);
trace = raw.data(:, 1, 1, ceil(params.sar.num_positions / 2));
plot(raw.fast_time * 1e6, real(trace), 'LineWidth', 1);
grid on;
xlabel('Fast time (us)');
ylabel('Real amplitude');
title('Example dechirped trace');

subplot(2, 2, 2);
range_profile = abs(range_data.data(:, 1, 1, ceil(params.sar.num_positions / 2)));
plot(range_data.range_axis, 20 * log10(range_profile ./ max(range_profile) + eps), 'LineWidth', 1);
grid on;
xlabel('Equivalent range (m)');
ylabel('Magnitude (dB)');
title('Example range profile');
ylim([-80, 5]);

subplot(2, 2, [3, 4]);
imagesc(image.x_axis, image.y_axis, image.magnitude_db);
set(gca, 'YDir', 'normal');
axis image;
if exist('turbo', 'file') || exist('turbo', 'builtin')
    colormap(turbo);
else
    colormap(jet);
end
colorbar;
caxis([-45, 0]);
xlabel('Cross-range x (m)');
ylabel('Range y (m)');
title('Focused near-field MIMO SAR image');
hold on;

for idx = 1:numel(targets)
    target_pos = targets(idx).position;
    plot(target_pos(1), target_pos(2), 'wo', 'MarkerSize', 8, 'LineWidth', 1.5);
end

plot(image.peak_position(1), image.peak_position(2), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
end
