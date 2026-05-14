# Near-field MIMO SAR Imaging Model

This repository contains a MATLAB simulation for near-field FMCW MIMO SAR imaging.
The first version models a one-dimensional linear SAR scan and reconstructs a
two-dimensional range/cross-range image with exact bistatic TX/RX path lengths.

## Run

Open MATLAB in this directory and run:

```matlab
main_sim_mimo_sar
```

The script will:

1. Load default radar, array, trajectory, and image-grid parameters.
2. Simulate dechirped FMCW data for multiple TX/RX channels over the SAR path.
3. Apply range FFT compression.
4. Reconstruct the near-field image with backprojection.
5. Plot an example trace, range profile, and focused image.

## Smoke test

Run this in MATLAB to verify that a noiseless single point target focuses within
the expected image/range-bin tolerance:

```matlab
tests/run_smoke_test
```

Run the PGA compensation integration test to verify that a known
RX/position-dependent MIMO phase error is estimated and compensated before
backprojection:

```matlab
tests/run_pga_compensation_test
```

For visual diagnostics, run:

```matlab
run_pga_mimo_compensation_experiment
```

By default this uses the configured point targets. To replace them with a
picture-derived sparse reflectivity scene, pass an image path:

```matlab
run_pga_mimo_compensation_experiment(true, '/path/to/target.png')
```

The image mode keeps about 3000 strongest dark pixels as real positive
scatterers, maps them to a 60 mm wide target in the imaging plane, and writes an
`image_target_reference.png` next to the experiment figures. The experiment then
compares clean, phase-degraded, and compensated images for a position-dependent
aperture phase-error case and a displacement-equivalent phase-error case. The
displacement case uses the first-order phase error implied by a small platform
perturbation; full unknown-trajectory motion compensation is a separate problem.

## Important parameters

Edit `config/default_params.m` to change:

- Radar waveform: center frequency, bandwidth, chirp time, ADC sample rate, samples per chirp.
- MIMO geometry: `params.array.tx_pos` and `params.array.rx_pos`.
- SAR path: `params.sar.platform_pos` or the start/stop/position count settings.
- Image grid: `params.image.x_axis`, `params.image.y_axis`, and fixed plane `params.image.z0`.
- Targets: `params.scene.targets`, where each row is `[x, y, z, complex_reflectivity]`.
- Noise and processing options: SNR, range window, zero-padding, and spreading-loss compensation.

## Model notes

For every platform position, TX, RX, and image pixel, the focusing model uses:

```matlab
L = norm(pixel - tx_pos) + norm(pixel - rx_pos);
R_equiv = L / 2;
```

The range-compressed signal is interpolated at `R_equiv`, then phase corrected
with the exact carrier phase associated with `L`. This avoids a far-field virtual
array approximation and is suitable for the first near-field modeling baseline.
