clear; clc;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
projectRoot = fileparts(repoRoot);
addpath(projectRoot);
addpath(genpath(repoRoot));

tmpImage = [tempname, '.png'];
img = ones(32, 32);
img(9:24, 13:20) = 0;
img(13:20, 8:25) = 0;
imwrite(img, tmpImage);

params = default_params();
params.radar.num_samples = 64;
params.processing.zero_pad_factor = 1;
params.processing.noise_snr_db = Inf;
params.processing.apply_spreading_loss = true;
params.processing.aperture_window = 'hann';
params.sar.platform_pos = [linspace(-0.002, 0.002, 5).', zeros(5, 1), zeros(5, 1)];
params.sar.num_positions = size(params.sar.platform_pos, 1);
params.scene.target_mode = 'image';
params.scene.image_path = tmpImage;
params.scene.image_num_scatterers = 30;
params.scene.image_target_width_m = 0.020;
params.scene.image_target_center = [0.0, 0.340];
params.image.x_axis = linspace(-0.020, 0.020, 21);
params.image.y_axis = linspace(0.320, 0.360, 21);
params.image.z0 = params.scene.z0;

targets = make_point_targets(params);
assert(numel(targets) == params.scene.image_num_scatterers, ...
    'Image target did not keep the requested number of scatterers.');

positions = reshape([targets.position], 3, []).';
reflectivity = [targets.reflectivity].';
assert(all(isfinite(positions(:))), 'Image target positions contain NaN or Inf.');
assert(all(isfinite(reflectivity(:))), 'Image target reflectivity contains NaN or Inf.');
assert(all(abs(imag(reflectivity)) < eps), 'Image target reflectivity should be real.');
assert(all(real(reflectivity) > 0), 'Image target reflectivity should be positive.');
assert(max(positions(:, 1)) - min(positions(:, 1)) <= params.scene.image_target_width_m + eps, ...
    'Image target x extent exceeds configured width.');

raw = simulate_fmcw_mimo_sar(params, targets);
rangeData = range_fft(raw, params);
image = backprojection_nearfield_mimo(rangeData, params);
assert(all(isfinite(image.magnitude(:))), 'Image target clean BP contains NaN or Inf.');
assert(image.peak_position(2) >= min(params.image.y_axis) && ...
    image.peak_position(2) <= max(params.image.y_axis), ...
    'Image target peak is outside the y image axis.');

delete(tmpImage);
fprintf('Image target test passed with %d scatterers.\n', numel(targets));
