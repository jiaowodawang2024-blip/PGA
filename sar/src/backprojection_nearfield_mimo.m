function image = backprojection_nearfield_mimo(range_data, params)
%BACKPROJECTION_NEARFIELD_MIMO Focus near-field MIMO SAR data with exact paths.

c = params.constants.c;
fc = params.radar.fc;

x_axis = params.image.x_axis;
y_axis = params.image.y_axis;
z0 = params.image.z0;
[x_grid, y_grid] = meshgrid(x_axis, y_axis);
z_grid = z0 * ones(size(x_grid));

platform_pos = params.sar.platform_pos;
tx_rel = params.array.tx_pos;
rx_rel = params.array.rx_pos;

num_positions = size(platform_pos, 1);
num_tx = size(tx_rel, 1);
num_rx = size(rx_rel, 1);
pos_window = aperture_weights(params, num_positions);
tx_window = make_imaging_window(get_processing_option(params, 'tx_window', 'rect'), num_tx);
rx_window = make_imaging_window(get_processing_option(params, 'rx_window', 'rect'), num_rx);
weight_sum = sum(pos_window) * sum(tx_window) * sum(rx_window);

focused = complex(zeros(size(x_grid)));

for pos_idx = 1:num_positions
    platform = platform_pos(pos_idx, :);
    w_pos = pos_window(pos_idx);

    for tx_idx = 1:num_tx
        tx_pos = platform + tx_rel(tx_idx, :);
        w_tx = tx_window(tx_idx);
        path_tx = sqrt( ...
            (x_grid - tx_pos(1)).^2 + ...
            (y_grid - tx_pos(2)).^2 + ...
            (z_grid - tx_pos(3)).^2);

        for rx_idx = 1:num_rx
            rx_pos = platform + rx_rel(rx_idx, :);
            w_rx = rx_window(rx_idx);
            path_rx = sqrt( ...
                (x_grid - rx_pos(1)).^2 + ...
                (y_grid - rx_pos(2)).^2 + ...
                (z_grid - rx_pos(3)).^2);

            path_length = path_tx + path_rx;
            equivalent_range = path_length / 2;

            sample = interp1( ...
                range_data.range_axis, ...
                range_data.data(:, rx_idx, tx_idx, pos_idx), ...
                equivalent_range, ...
                params.processing.interpolation, ...
                0);

            phase_correction = exp(1j * 2 * pi * fc * path_length / c);

            if params.processing.apply_spreading_loss
                sample = sample .* max(path_tx .* path_rx, eps);
            end

            focused = focused + w_pos * w_tx * w_rx * sample .* phase_correction;
        end
    end
end

if weight_sum > eps
    focused = focused * (num_positions * num_tx * num_rx / weight_sum);
end

magnitude = abs(focused);
[~, peak_idx] = max(magnitude(:));
[peak_y_idx, peak_x_idx] = ind2sub(size(magnitude), peak_idx);

image.complex = focused;
image.magnitude = magnitude;
image.magnitude_db = 20 * log10(magnitude ./ max(magnitude(:)) + eps);
image.x_axis = x_axis;
image.y_axis = y_axis;
image.z0 = z0;
image.peak_position = [x_axis(peak_x_idx), y_axis(peak_y_idx), z0];
end

function value = get_processing_option(params, name, defaultValue)
value = defaultValue;
if isfield(params, 'processing') && isfield(params.processing, name) && ~isempty(params.processing.(name))
    value = params.processing.(name);
end
end

function weights = aperture_weights(params, num_positions)
if isfield(params, 'sar') && isfield(params.sar, 'platform_weights') && ...
        ~isempty(params.sar.platform_weights)
    weights = params.sar.platform_weights(:);
    if numel(weights) ~= num_positions
        error('backprojection_nearfield_mimo:BadPlatformWeights', ...
            'params.sar.platform_weights must have one value per platform position.');
    end
    return;
end

weights = make_imaging_window(get_processing_option(params, 'aperture_window', 'rect'), num_positions);
end

function window = make_imaging_window(name, n)
if n <= 1
    window = ones(n, 1);
    return;
end

switch lower(char(name))
    case {'hann', 'hanning'}
        k = (0:n - 1).';
        window = 0.5 - 0.5 * cos(2 * pi * k / (n - 1));
    case {'rect', 'rectangular', 'none'}
        window = ones(n, 1);
    otherwise
        error('backprojection_nearfield_mimo:BadWindow', ...
            'Unsupported imaging window: %s', char(name));
end
end
