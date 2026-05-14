clear; clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));

params = default_params();
params.processing.noise_snr_db = Inf;
params.scene.targets = [0.020, 0.330, 0, 1.0];
params.image.x_axis = linspace(-0.04, 0.04, 81);
params.image.y_axis = linspace(0.28, 0.38, 101);

targets = make_point_targets(params);
raw = simulate_fmcw_mimo_sar(params, targets);
range_data = range_fft(raw, params);
image = backprojection_nearfield_mimo(range_data, params);

truth = targets(1).position;
dx = abs(image.peak_position(1) - truth(1));
dy = abs(image.peak_position(2) - truth(2));
grid_dx = mean(diff(params.image.x_axis));
grid_dy = mean(diff(params.image.y_axis));
range_bin = mean(diff(range_data.range_axis));

assert(dx <= grid_dx, 'Peak x error %.6f exceeds one grid cell %.6f.', dx, grid_dx);
assert(dy <= max(grid_dy, range_bin), ...
    'Peak y error %.6f exceeds tolerance %.6f.', dy, max(grid_dy, range_bin));

fprintf('Smoke test passed. Peak: x = %.4f m, y = %.4f m\n', image.peak_position(1), image.peak_position(2));
