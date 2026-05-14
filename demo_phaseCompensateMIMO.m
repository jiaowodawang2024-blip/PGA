% Demo for THz MIMO RX phase compensation during 2-D motion scanning.
%
% The synthetic data below has layout data(rx, tx, freq, y, x). It injects a
% known RX-dependent, position-dependent phase jitter and then estimates it
% using only channel coherence against a reference RX channel.

clear; clc; close all;
rng(7);

nRx = 6;
nTx = 4;
nFreq = 96;
nY = 48;
nX = 64;
refRx = 1;

[yy, xx] = ndgrid(linspace(-1, 1, nY), linspace(-1, 1, nX));
scenePhase = 2.0 * xx + 1.2 * yy + 0.7 * exp(-6 * (xx.^2 + yy.^2));
sceneAmp = 1.0 + 0.45 * exp(-18 * ((xx - 0.25).^2 + (yy + 0.18).^2));

freqPhase = reshape(linspace(0, 8 * pi, nFreq), [1 1 nFreq 1 1]);
txPhase = reshape(linspace(0, 0.6 * pi, nTx), [1 nTx 1 1 1]);
baseScene = reshape(sceneAmp .* exp(1i * scenePhase), [1 1 1 nY nX]);
dataClean = baseScene .* exp(1i * freqPhase) .* exp(1i * txPhase);
dataClean = repmat(dataClean, [nRx 1 1 1 1]);

rxStaticPhase = zeros(nRx, 1, 1);
phaseTrue = zeros(nRx, nY, nX);
for rx = 1:nRx
    phaseTrue(rx, :, :) = rxStaticPhase(rx) ...
        + 0.45 * sin(2 * pi * (0.8 * xx + 0.15 * rx)) ...
        + 0.32 * cos(2 * pi * (0.6 * yy - 0.08 * rx)) ...
        + 0.18 * sin(2 * pi * (xx + yy));
end
phaseTrue(refRx, :, :) = 0;

jitter = exp(1i * reshape(phaseTrue, [nRx 1 1 nY nX]));
noise = 0.12 / sqrt(2) * (randn(size(dataClean)) + 1i * randn(size(dataClean)));
dataJittered = dataClean .* jitter + noise;

opts = struct();
opts.refRx = refRx;
opts.dimOrder = 'rx_tx_freq_y_x';
opts.method = 'pga';
opts.freqIdx = 12:84;
opts.pgaNumRangeBins = 8;
opts.smoothWindow = [7 7];
opts.unwrapSpatial = true;
opts.removeMeanPhase = true;
opts.minMagnitudePercentile = 10;

[dataCorr, phaseEst, info] = phaseCompensateMIMO(dataJittered, opts);

phaseEstStd = squeeze(phaseEst(:, 1, 1, :, :));
phaseErrBefore = angle(exp(1i * (phaseTrue - phaseTrue(refRx, :, :))));
phaseErrAfter = angle(exp(1i * (phaseEstStd - phaseErrBefore)));
rmseRad = sqrt(mean(phaseErrAfter(:).^2));

imgBeforeComplex = squeeze(sum(sum(sum(dataJittered(:, :, opts.freqIdx, :, :), 1), 2), 3));
imgAfterComplex = squeeze(sum(sum(sum(dataCorr(:, :, opts.freqIdx, :, :), 1), 2), 3));
imgBefore = abs(imgBeforeComplex);
imgAfter = abs(imgAfterComplex);
imgBeforeDb = normalizeImageDb(imgBefore);
imgAfterDb = normalizeImageDb(imgAfter);
imgGainDb = imgAfterDb - imgBeforeDb;
contrastBefore = std(imgBefore(:)) / mean(imgBefore(:));
contrastAfter = std(imgAfter(:)) / mean(imgAfter(:));

fprintf('Mean coherence before: %.4f\n', mean(info.coherenceBefore, 'omitnan'));
fprintf('Mean coherence after : %.4f\n', mean(info.coherenceAfter, 'omitnan'));
fprintf('Phase estimate RMSE  : %.4f rad\n', rmseRad);
fprintf('Image contrast before: %.4f\n', contrastBefore);
fprintf('Image contrast after : %.4f\n', contrastAfter);

rxShow = min(3, nRx);
figure('Name', 'MIMO RX Phase Compensation Demo', 'Color', 'w');

subplot(2, 3, 1);
imagesc(squeeze(phaseErrBefore(rxShow, :, :)));
axis image; colorbar;
xlabel('Scan x index');
ylabel('Scan y index');
title(sprintf('True phase, RX %d', rxShow));

subplot(2, 3, 2);
imagesc(squeeze(phaseEstStd(rxShow, :, :)));
axis image; colorbar;
xlabel('Scan x index');
ylabel('Scan y index');
title(sprintf('Estimated phase, RX %d', rxShow));

subplot(2, 3, 3);
imagesc(squeeze(phaseErrAfter(rxShow, :, :)));
axis image; colorbar;
xlabel('Scan x index');
ylabel('Scan y index');
title('Residual phase error');

subplot(2, 3, 4);
bar([info.coherenceBefore(:), info.coherenceAfter(:)]);
grid on;
xlabel('RX channel');
ylabel('Coherence to reference');
legend('Before', 'After', 'Location', 'best');
title('Channel coherence');

subplot(2, 3, 5);
imagesc(imgBefore);
axis image; colorbar;
xlabel('Scan x index');
ylabel('Scan y index');
title('Coherent sum before');

subplot(2, 3, 6);
imagesc(imgAfter);
axis image; colorbar;
xlabel('Scan x index');
ylabel('Scan y index');
title('Coherent sum after');

figure('Name', 'MIMO Image Before and After Phase Compensation', 'Color', 'w');

subplot(1, 3, 1);
imagesc(imgBeforeDb);
axis image; colorbar;
caxis([-40 0]);
xlabel('Scan x index');
ylabel('Scan y index');
title('Image before compensation (dB)');

subplot(1, 3, 2);
imagesc(imgAfterDb);
axis image; colorbar;
caxis([-40 0]);
xlabel('Scan x index');
ylabel('Scan y index');
title('Image after compensation (dB)');

subplot(1, 3, 3);
imagesc(imgGainDb);
axis image; colorbar;
xlabel('Scan x index');
ylabel('Scan y index');
title('Image gain after compensation (dB)');

function imgDb = normalizeImageDb(img)
img = abs(img);
img = img ./ max(img(:) + eps);
imgDb = 20 * log10(img + eps);
end
