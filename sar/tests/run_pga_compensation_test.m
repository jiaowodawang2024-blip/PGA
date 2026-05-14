clear; clc;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
projectRoot = fileparts(repoRoot);
addpath(projectRoot);
addpath(genpath(repoRoot));

results = run_pga_mimo_compensation_experiment(false);

assert(isfield(results, 'imageClean'), 'Missing clean image.');
assert(isfield(results.positionRxPhase, 'imageDegraded'), 'Missing degraded image.');
assert(isfield(results.positionRxPhase, 'imageCompensated'), 'Missing compensated image.');
assert(all(isfinite(results.cleanData(:))), 'cleanData contains NaN or Inf.');
assert(all(isfinite(results.positionRxPhase.degradedData(:))), 'degradedData contains NaN or Inf.');
assert(all(isfinite(results.positionRxPhase.compensatedData(:))), 'compensatedData contains NaN or Inf.');
assert(results.params.sar.num_x_positions == 41, 'Unexpected 2-D aperture x size.');
assert(results.params.sar.num_y_positions == 41, 'Unexpected 2-D aperture y size.');
assert(results.params.sar.num_positions == 41 * 41, 'Unexpected 2-D aperture position count.');

lambda = results.params.radar.lambda;
platformPos = results.params.sar.platform_pos;
xUnique = unique(round(platformPos(:, 1), 12));
yUnique = unique(round(platformPos(:, 2), 12));
assert(abs(median(diff(xUnique)) - lambda / 2) < 1e-12, ...
    'Aperture x spacing is not lambda/2.');
assert(abs(median(diff(yUnique)) - lambda / 2) < 1e-12, ...
    'Aperture y spacing is not lambda/2.');

txPos = results.params.array.tx_pos;
rxPos = results.params.array.rx_pos;
phaseCenters = zeros(size(txPos, 1) * size(rxPos, 1), 3);
idx = 1;
for txIdx = 1:size(txPos, 1)
    for rxIdx = 1:size(rxPos, 1)
        phaseCenters(idx, :) = 0.5 * (txPos(txIdx, :) + rxPos(rxIdx, :));
        idx = idx + 1;
    end
end
centerX = unique(round(phaseCenters(:, 1), 12));
centerY = unique(round(phaseCenters(:, 2), 12));
assert(numel(centerX) == 4 && numel(centerY) == 4, ...
    'Virtual phase centers are not a 4 x 4 grid.');
assert(abs(median(diff(centerX)) - lambda / 2) < 1e-12, ...
    'Virtual phase-center x spacing is not lambda/2.');
assert(abs(median(diff(centerY)) - lambda / 2) < 1e-12, ...
    'Virtual phase-center y spacing is not lambda/2.');

clean = results.metrics.clean;
posDeg = results.positionRxPhase.metrics.degraded;
posComp = results.positionRxPhase.metrics.compensated;
motDeg = results.motionEquivalent.metrics.degraded;
motComp = results.motionEquivalent.metrics.compensated;

assert(clean.peakValue > 0, 'Clean image peak is zero.');
assert(posDeg.imageEntropy >= clean.imageEntropy || posDeg.ISLR >= clean.ISLR || ...
    posDeg.energyConcentration <= clean.energyConcentration, ...
    'Degraded image did not show a measurable focus degradation.');

positionImproved = posComp.PSLR < posDeg.PSLR || ...
    posComp.ISLR < posDeg.ISLR || ...
    posComp.imageEntropy < posDeg.imageEntropy || ...
    posComp.energyConcentration > posDeg.energyConcentration;
motionImproved = motComp.ISLR < motDeg.ISLR && ...
    motComp.imageEntropy < motDeg.imageEntropy && ...
    motComp.peakValue > 2.0 * motDeg.peakValue;

if ~(positionImproved && motionImproved)
    debugFile = fullfile(results.outputDir, 'debug_test_failure.mat');
    save(debugFile, 'results');
    warning('run_pga_compensation_test:NoImprovement', ...
        'Metrics did not improve in all scenarios. Debug results saved to %s', debugFile);
end

assert(positionImproved, 'Position/RX phase compensation did not improve degraded metrics.');
assert(motionImproved, 'Motion-equivalent phase compensation did not improve degraded metrics.');

fprintf('PGA MIMO compensation test passed.\n');
