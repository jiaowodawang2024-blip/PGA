clear; clc; close all;

addpath(genpath(fileparts(mfilename('fullpath'))));

params = default_params();

targets = make_point_targets(params);
raw = simulate_fmcw_mimo_sar(params, targets);
range_data = range_fft(raw, params);
image = backprojection_nearfield_mimo(range_data, params);

plot_results(raw, range_data, image, targets, params);

fprintf('MIMO SAR simulation complete.\n');
fprintf('Image peak: x = %.4f m, y = %.4f m\n', image.peak_position(1), image.peak_position(2));
