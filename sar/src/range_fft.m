function range_data = range_fft(raw, params)
%RANGE_FFT Compress FMCW fast-time data into equivalent monostatic range bins.

num_samples = params.radar.num_samples;
nfft = 2^nextpow2(num_samples * params.processing.zero_pad_factor);

window = make_window(params.processing.range_window, num_samples);
windowed = raw.data .* reshape(window, [], 1, 1, 1);
spectrum = fft(windowed, nfft, 1);

positive_bins = 1:floor(nfft / 2) + 1;
spectrum = spectrum(positive_bins, :, :, :);

frequency_axis = (positive_bins - 1).' * params.radar.sample_rate / nfft;
range_axis = params.constants.c * frequency_axis / (2 * params.radar.slope);

range_data.data = spectrum;
range_data.range_axis = range_axis;
range_data.frequency_axis = frequency_axis;
range_data.nfft = nfft;
end

function window = make_window(name, num_samples)
switch lower(name)
    case {'hann', 'hanning'}
        if num_samples == 1
            window = 1;
        else
            n = (0:num_samples - 1).';
            window = 0.5 - 0.5 * cos(2 * pi * n / (num_samples - 1));
        end
    case {'rect', 'rectangular', 'none'}
        window = ones(num_samples, 1);
    otherwise
        error('Unsupported range window: %s', name);
end
end
