function raw = simulate_fmcw_mimo_sar(params, targets)
%SIMULATE_FMCW_MIMO_SAR Simulate dechirped FMCW data for a MIMO SAR scan.

rng(params.processing.random_seed);

c = params.constants.c;
fc = params.radar.fc;
slope = params.radar.slope;
fs = params.radar.sample_rate;
num_samples = params.radar.num_samples;

platform_pos = params.sar.platform_pos;
tx_rel = params.array.tx_pos;
rx_rel = params.array.rx_pos;

num_positions = size(platform_pos, 1);
num_tx = size(tx_rel, 1);
num_rx = size(rx_rel, 1);
t_fast = (0:num_samples - 1).' / fs;

data = complex(zeros(num_samples, num_rx, num_tx, num_positions));

for pos_idx = 1:num_positions
    platform = platform_pos(pos_idx, :);

    for tx_idx = 1:num_tx
        tx_pos = platform + tx_rel(tx_idx, :);

        for rx_idx = 1:num_rx
            rx_pos = platform + rx_rel(rx_idx, :);
            trace = complex(zeros(num_samples, 1));

            for target_idx = 1:numel(targets)
                target_pos = targets(target_idx).position;
                sigma = targets(target_idx).reflectivity;

                path_tx = norm(target_pos - tx_pos);
                path_rx = norm(target_pos - rx_pos);
                path_length = path_tx + path_rx;
                tau = path_length / c;
                beat_frequency = slope * tau;

                phase0 = -2 * pi * fc * tau + pi * slope * tau^2;
                amplitude = sigma;

                if params.processing.apply_spreading_loss
                    amplitude = amplitude / max(path_tx * path_rx, eps);
                end

                trace = trace + amplitude .* exp(1j * (2 * pi * beat_frequency .* t_fast + phase0));
            end

            data(:, rx_idx, tx_idx, pos_idx) = trace;
        end
    end
end

if isfinite(params.processing.noise_snr_db)
    signal_power = mean(abs(data(:)).^2);
    noise_power = signal_power / (10^(params.processing.noise_snr_db / 10));
    noise = sqrt(noise_power / 2) .* (randn(size(data)) + 1j * randn(size(data)));
    data = data + noise;
end

raw.data = data;
raw.fast_time = t_fast;
raw.dimensions = struct( ...
    'num_samples', num_samples, ...
    'num_rx', num_rx, ...
    'num_tx', num_tx, ...
    'num_positions', num_positions);
end
