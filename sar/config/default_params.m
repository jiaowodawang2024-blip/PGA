function params = default_params()
%DEFAULT_PARAMS Default parameters for near-field FMCW MIMO SAR simulation.

c = 299792458;

params.constants.c = c;

params.radar.fc = 300e9;
params.radar.bandwidth = 20e9;
params.radar.chirp_time = 100e-6;
params.radar.sample_rate = 5e6;
params.radar.num_samples = 512;
params.radar.slope = params.radar.bandwidth / params.radar.chirp_time;
params.radar.lambda = c / params.radar.fc;

lambda = params.radar.lambda;

% Antenna coordinates are relative to the SAR platform reference point.
params.array.tx_pos = [ ...
    -0.75 * lambda, 0, 0; ...
     0.75 * lambda, 0, 0];
params.array.rx_pos = [ ...
    -0.75 * lambda, 0, 0; ...
    -0.25 * lambda, 0, 0; ...
     0.25 * lambda, 0, 0; ...
     0.75 * lambda, 0, 0];

params.sar.num_positions = 81;
params.sar.x_start = -0.08;
params.sar.x_stop = 0.08;
params.sar.y = 0;
params.sar.z = 0;

x_positions = linspace(params.sar.x_start, params.sar.x_stop, params.sar.num_positions).';
params.sar.platform_pos = [ ...
    x_positions, ...
    params.sar.y * ones(params.sar.num_positions, 1), ...
    params.sar.z * ones(params.sar.num_positions, 1)];

params.scene.z0 = 0;
params.scene.targets = [ ...
    -0.025, 0.300, 0, 1.00 * exp(1j * 0.0); ...
     0.030, 0.355, 0, 0.80 * exp(1j * 0.7); ...
     0.000, 0.430, 0, 0.55 * exp(1j * 1.4)];

params.image.x_axis = linspace(-0.10, 0.10, 161);
params.image.y_axis = linspace(0.20, 0.50, 181);
params.image.z0 = params.scene.z0;

params.processing.range_window = 'hann';
params.processing.zero_pad_factor = 4;
params.processing.apply_spreading_loss = true;
params.processing.noise_snr_db = 35;
params.processing.random_seed = 7;
params.processing.interpolation = 'linear';
end
