function results = run_pga_mimo_compensation_experiment(makePlots, imagePath)
%RUN_PGA_MIMO_COMPENSATION_EXPERIMENT MIMO RX phase-error autofocus experiment.
%
%   results = run_pga_mimo_compensation_experiment()
%   results = run_pga_mimo_compensation_experiment(false)
%   results = run_pga_mimo_compensation_experiment(true, imagePath)
%
% Engineering data dimensions used by this project:
%   rawClean.data(scanFreqSample, rx, tx, scanPos)
%   rangeClean.data(rangeBin, rx, tx, scanPos)
%
% The user-facing preferred notation rawData(scanPos, tx, rx, freqSample)
% is converted here by explicit permute/reshape calls. Phase errors are
% injected into complex range data, not into focused images. All phases are
% stored in radians; plots convert to degrees only for display.

if nargin < 1 || isempty(makePlots)
    makePlots = true;
end
if nargin < 2
    imagePath = '';
end

thisDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(thisDir);
addpath(repoRoot);
addpath(genpath(thisDir));

outDir = fullfile(thisDir, 'results', 'pga_mimo_compensation');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

params = experimentParams(imagePath);
targets = make_point_targets(params);
targetInfo = targetMetadata(params, targets);

rawClean = simulate_fmcw_mimo_sar(params, targets);
rangeClean = range_fft(rawClean, params);
imageClean = backprojection_nearfield_mimo(rangeClean, params);

positionScenario = runPositionRxScenario(rangeClean, imageClean, targets, params);
motionScenario = runMotionEquivalentScenario(rangeClean, imageClean, targets, params);

results = struct();
results.params = params;
results.targets = targets;
results.targetInfo = targetInfo;
results.cleanData = rangeClean.data;
results.imageClean = imageClean;
results.metrics.clean = imageQualityMetrics(imageClean);
results.positionRxPhase = positionScenario;
results.motionEquivalent = motionScenario;
results.outputDir = outDir;
results.psfDiagnostics = psfDiagnostics(params, imageClean);

printMetricsTable('Common aperture + fixed RX phase error', ...
    results.metrics.clean, positionScenario.metrics.degraded, positionScenario.metrics.compensated);
printConvergenceTable('Common aperture + fixed RX phase error convergence', positionScenario.iterationMetrics);
printMetricsTable('Motion-equivalent phase error', ...
    results.metrics.clean, motionScenario.metrics.degraded, motionScenario.metrics.compensated);
printConvergenceTable('Motion-equivalent phase error convergence', motionScenario.iterationMetrics);

if ~isImproved(positionScenario.metrics.degraded, positionScenario.metrics.compensated)
    warning('run_pga_mimo_compensation_experiment:NoPositionImprovement', ...
        'Compensated position/RX metrics did not improve over degraded metrics.');
end
if ~isImproved(motionScenario.metrics.degraded, motionScenario.metrics.compensated)
    warning('run_pga_mimo_compensation_experiment:NoMotionImprovement', ...
        'Compensated motion metrics did not improve over degraded metrics.');
end

save(fullfile(outDir, 'pga_mimo_compensation_results.mat'), 'results');

if makePlots
    if strcmpi(targetInfo.mode, 'image')
        plotImageTargetReference(outDir, targets, targetInfo);
    end
    plotScenario(outDir, imageClean, positionScenario, ...
        'Common aperture + fixed RX phase error', 'position_rx_phase');
    plotScenario(outDir, imageClean, motionScenario, ...
        'Motion-equivalent phase error', 'motion_equivalent');
    plotConvergence(outDir, positionScenario, ...
        'Common aperture + fixed RX phase error convergence', 'position_rx_phase_convergence');
    plotConvergence(outDir, motionScenario, ...
        'Motion-equivalent phase error convergence', 'motion_equivalent_convergence');
end
end

function params = experimentParams(imagePath)
if nargin < 1
    imagePath = '';
end

params = default_params();
params.radar.num_samples = 256;
params.radar.bandwidth = 40e9;
params.radar.slope = params.radar.bandwidth / params.radar.chirp_time;
params.processing.zero_pad_factor = 2;
params.processing.noise_snr_db = Inf;
params.processing.random_seed = 19;
params.processing.apply_spreading_loss = true;
params.processing.aperture_window = 'hann';
params.processing.tx_window = 'rect';
params.processing.rx_window = 'rect';

lambda = params.radar.lambda;
d = lambda / 2;

[txGridX, txGridY] = meshgrid([-1 1] * lambda, [-1 1] * lambda);
[rxGridX, rxGridY] = meshgrid([-0.5 0.5] * lambda, [-0.5 0.5] * lambda);
params.array.tx_pos = [txGridX(:), txGridY(:), zeros(numel(txGridX), 1)];
params.array.rx_pos = [rxGridX(:), rxGridY(:), zeros(numel(rxGridX), 1)];

scanNx = 41;
scanNy = 41;
x = ((0:scanNx - 1) - (scanNx - 1) / 2) * d;
y = ((0:scanNy - 1) - (scanNy - 1) / 2) * d;
[X, Y] = meshgrid(x, y);
params.sar.num_x_positions = scanNx;
params.sar.num_y_positions = scanNy;
params.sar.aperture_step = d;
params.sar.num_positions = numel(X);
params.sar.platform_pos = [X(:), Y(:), zeros(numel(X), 1)];
params.sar.platform_weights = makeSeparableHannWeights(scanNy, scanNx);

params.scene.targets = [ ...
    -0.015, 0.335, 0, 1.00; ...
     0.018, 0.350, 0, 0.65 * exp(1j * 0.6); ...
     0.036, 0.322, 0, 0.35 * exp(1j * 1.2)];
if ~isempty(imagePath)
    params.scene.target_mode = 'image';
    params.scene.image_path = char(imagePath);
    params.scene.image_num_scatterers = 3000;
    params.scene.image_target_width_m = 0.060;
else
    params.scene.target_mode = 'points';
end

params.image.x_axis = linspace(-0.060, 0.060, 101);
params.image.y_axis = linspace(0.290, 0.390, 101);
params.image.z0 = params.scene.z0;
if strcmpi(params.scene.target_mode, 'image')
    params.scene.image_target_center = [mean(params.image.x_axis), mean(params.image.y_axis)];
end

params.motion_error.amplitude_mm = 0.45;
params.motion_error.period_positions = 10;
params.motion_error.num_pga_iterations = 3;
end

function scenario = runPositionRxScenario(rangeClean, imageClean, targets, params)
nRx = size(rangeClean.data, 2);
nPos = size(rangeClean.data, 4);
refRx = 1;

phiCommonTrue = makeCommonAperturePhase(nPos, 1.15);
phiRxTrue = makeFixedRxPhase(nRx, refRx, 0.75);
phiTrue = bsxfun(@plus, phiCommonTrue(:).', phiRxTrue(:));

rangeDegraded = applyPhaseMap(rangeClean, phiTrue);
imageDegraded = backprojection_nearfield_mimo(rangeDegraded, params);

opts = compensationOptions(params, refRx, 'rx_phase');
[opts.nominalGradientX, opts.nominalGradientY] = ...
    makeNominalGradients(params, dominantTargetPosition(targets));
[rangeCompensated, phiHat, info, iterationImages, iterationMetrics] = compensateRangeData(rangeDegraded, opts, params);
imageCompensated = backprojection_nearfield_mimo(rangeCompensated, params);

scenario = buildScenario(rangeClean, rangeDegraded, rangeCompensated, ...
    imageClean, imageDegraded, imageCompensated, phiTrue, phiCommonTrue, ...
    phiRxTrue, phiHat, info, targets);
scenario.iterationImages = iterationImages;
scenario.iterationMetrics = iterationMetrics;
end

function scenario = runMotionEquivalentScenario(rangeClean, imageClean, targets, params)
nRx = size(rangeClean.data, 2);
nPos = size(rangeClean.data, 4);
refRx = 1;

[dx, dy] = makePlatformDisplacement(params, nPos);
phiMotionTrue = makeMotionEquivalentPhase(params, targets, dx, dy);
phiTrue = repmat(phiMotionTrue(:).', nRx, 1);

rangeDegraded = applyPhaseMap(rangeClean, phiTrue);
imageDegraded = backprojection_nearfield_mimo(rangeDegraded, params);

opts = compensationOptions(params, refRx, 'pga');
opts.phaseModel = 'aperture_common';
[opts.nominalGradientX, opts.nominalGradientY] = ...
    makeNominalGradients(params, dominantTargetPosition(targets));
[rangeBaselinePga, phiBaselinePga, infoBaselinePga, baselineImages, baselineMetrics] = ...
    compensateRangeData(rangeDegraded, opts, params);

optsGeom = geometryCompensationOptions(params);
[rangeCompensated, phiHat, info, iterationImages, iterationMetrics] = ...
    compensateMotionGeometry(rangeDegraded, optsGeom, params);
imageCompensated = backprojection_nearfield_mimo(rangeCompensated, params);

scenario = buildScenario(rangeClean, rangeDegraded, rangeCompensated, ...
    imageClean, imageDegraded, imageCompensated, phiTrue, phiMotionTrue, ...
    zeros(nRx, 1), phiHat, info, targets);
scenario.platformDisplacement = [dx(:), dy(:)];
scenario.iterationImages = iterationImages;
scenario.iterationMetrics = iterationMetrics;
scenario.baselinePga = struct( ...
    'compensatedData', rangeBaselinePga.data, ...
    'phiHat', phiBaselinePga, ...
    'info', infoBaselinePga, ...
    'iterationImages', baselineImages, ...
    'iterationMetrics', baselineMetrics, ...
    'imageCompensated', backprojection_nearfield_mimo(rangeBaselinePga, params), ...
    'metrics', imageQualityMetrics(backprojection_nearfield_mimo(rangeBaselinePga, params)));
end

function opts = compensationOptions(params, refRx, method)
opts = struct();
opts.method = method;
opts.pgaInputDomain = 'range';
opts.numIter = params.motion_error.num_pga_iterations;
opts.refRx = refRx;
opts.selectedRangeBins = [];
opts.selectedImageRegion = [];
opts.enableAmplitudeEq = true;
opts.unwrapPhase = true;
opts.dimOrder = 'rx_tx_freq_y_x';
opts.phaseModel = 'mimo_aperture';
opts.pgaNumRangeBins = 8;
opts.smoothWindow = [1 5];
opts.removeMeanPhase = true;
opts.removeLinearPhase = true;
opts.minMagnitudePercentile = 0;
end

function opts = geometryCompensationOptions(params)
opts = struct();
opts.numIter = params.motion_error.num_pga_iterations;
opts.peakSearchHalfWindow = [0.030, 0.030];
opts.phaseSmoothWindow = 1;
opts.unwrapPhase = true;
opts.minCoherenceFraction = 0.05;
opts.bestMetric = 'focused_quality';
end

function [rangeCompensated, phiHat, info, iterationImages, iterationMetrics] = ...
    compensateMotionGeometry(rangeData, opts, params)
rangeIter = rangeData;
nRx = size(rangeData.data, 2);
nPos = size(rangeData.data, 4);
phaseTotal = zeros(1, nPos);

iterationImages = repmat(backprojection_nearfield_mimo(rangeData, params), opts.numIter, 1);
iterationMetrics = repmat(imageQualityMetrics(iterationImages(1)), opts.numIter, 1);
iterationRanges = cell(opts.numIter, 1);
phaseSteps = zeros(opts.numIter, nPos);
focusPoints = zeros(opts.numIter, 3);
coherence = zeros(opts.numIter, nPos);

for iterIdx = 1:opts.numIter
    imageCurrent = backprojection_nearfield_mimo(rangeIter, params);
    focusPoint = selectFocusPoint(imageCurrent);
    [phaseStep, coh] = estimateGeometryCommonPhase(rangeIter, params, focusPoint, opts);
    rangeIter = applyCommonPhaseToRange(rangeIter, phaseStep);
    phaseTotal = phaseTotal + phaseStep(:).';

    iterationRanges{iterIdx} = rangeIter;
    iterationImages(iterIdx) = backprojection_nearfield_mimo(rangeIter, params);
    iterationMetrics(iterIdx) = imageQualityMetrics(iterationImages(iterIdx));
    phaseSteps(iterIdx, :) = phaseStep(:).';
    focusPoints(iterIdx, :) = focusPoint;
    coherence(iterIdx, :) = coh(:).';
end

bestIdx = selectBestIteration(iterationMetrics);
if bestIdx <= 0
    bestIdx = opts.numIter;
end
rangeCompensated = iterationRanges{bestIdx};
phiCommon = sum(phaseSteps(1:bestIdx, :), 1);
phiHat = reshape(repmat(phiCommon, [nRx 1]), [nRx 1 1 1 nPos]);

info = struct();
info.method = 'geometry_common_phase';
info.bestIteration = bestIdx;
info.bestIterationReason = 'minimum ISLR, then lower mainlobe width and entropy';
info.phaseSteps = phaseSteps;
info.phaseTotal = phiCommon;
info.phiCommonHat = phiCommon;
info.phiRxHat = zeros(nRx, 1);
info.focusPoints = focusPoints;
info.coherenceByPosition = coherence;
info.iterations = struct([]);
for iterIdx = 1:opts.numIter
    info.iterations(iterIdx).phaseStep = phaseSteps(iterIdx, :);
    info.iterations(iterIdx).phaseTotal = sum(phaseSteps(1:iterIdx, :), 1);
    info.iterations(iterIdx).focusPoint = focusPoints(iterIdx, :);
    info.iterations(iterIdx).coherenceByPosition = coherence(iterIdx, :);
end
end

function focusPoint = selectFocusPoint(image)
focusPoint = image.peak_position;
end

function [phaseCommon, coherence] = estimateGeometryCommonPhase(rangeData, params, focusPoint, opts)
c = params.constants.c;
fc = params.radar.fc;
platformPos = params.sar.platform_pos;
txRel = params.array.tx_pos;
rxRel = params.array.rx_pos;
nPos = size(rangeData.data, 4);
coherent = complex(zeros(1, nPos));

for posIdx = 1:nPos
    platform = platformPos(posIdx, :);
    acc = 0;
    for txIdx = 1:size(txRel, 1)
        txPos = platform + txRel(txIdx, :);
        pathTx = norm(focusPoint - txPos);
        for rxIdx = 1:size(rxRel, 1)
            rxPos = platform + rxRel(rxIdx, :);
            pathRx = norm(focusPoint - rxPos);
            pathLength = pathTx + pathRx;
            equivalentRange = pathLength / 2;
            sample = interp1(rangeData.range_axis, ...
                rangeData.data(:, rxIdx, txIdx, posIdx), ...
                equivalentRange, params.processing.interpolation, 0);
            if params.processing.apply_spreading_loss
                sample = sample .* max(pathTx .* pathRx, eps);
            end
            acc = acc + sample .* exp(1j * 2 * pi * fc * pathLength / c);
        end
    end
    coherent(posIdx) = acc;
end

coherence = abs(coherent);
phaseCommon = angle(coherent);
valid = coherence >= opts.minCoherenceFraction * max(coherence);
if nnz(valid) >= 2
    phaseCommon(~valid) = interp1(find(valid), phaseCommon(valid), find(~valid), 'linear', 'extrap');
end
if opts.unwrapPhase
    phaseCommon = unwrap(phaseCommon);
end
if opts.phaseSmoothWindow > 1
    kernel = ones(1, opts.phaseSmoothWindow) ./ opts.phaseSmoothWindow;
    phaseCommon = conv(phaseCommon, kernel, 'same');
end
phaseCommon = phaseCommon - mean(phaseCommon);
end

function rangeOut = applyCommonPhaseToRange(rangeIn, phaseCommon)
rangeOut = rangeIn;
rangeOut.data = bsxfun(@times, rangeIn.data, ...
    reshape(exp(-1i * phaseCommon(:).'), [1 1 1 numel(phaseCommon)]));
end

function [rangeCompensated, phiHat, info, iterationImages, iterationMetrics] = compensateRangeData(rangeData, opts, params)
% rangeData.data is [rangeBin, rx, tx, scanPos].
% phaseCompensateMIMO expects [rx, tx, freq/range, y, x].
nRange = size(rangeData.data, 1);
nRx = size(rangeData.data, 2);
nTx = size(rangeData.data, 3);
nPos = size(rangeData.data, 4);
[nY, nX] = apertureGridSize(params, nPos);

compInput = permute(rangeData.data, [2 3 1 4]);
compInput = reshape(compInput, [nRx nTx nRange nY nX]);
[compOutput, phiHat, info] = phaseCompensateMIMO(compInput, opts);

compOutput = reshape(compOutput, [nRx nTx nRange nPos]);
rangeCompensated = rangeData;
rangeCompensated.data = permute(compOutput, [3 1 2 4]);

iterationImages = [];
iterationMetrics = [];
if isfield(info, 'iterations') && ~isempty(info.iterations)
    iterationRanges = cell(numel(info.iterations), 1);
    iterationImages = repmat(backprojection_nearfield_mimo(rangeCompensated, params), numel(info.iterations), 1);
    iterationMetrics = repmat(imageQualityMetrics(iterationImages(1)), numel(info.iterations), 1);
    for iterIdx = 1:numel(info.iterations)
        if ~isfield(info.iterations(iterIdx), 'dataStd') || isempty(info.iterations(iterIdx).dataStd)
            continue;
        end
        iterOutput = reshape(info.iterations(iterIdx).dataStd, [nRx nTx nRange nPos]);
        iterRange = rangeData;
        iterRange.data = permute(iterOutput, [3 1 2 4]);
        iterationRanges{iterIdx} = iterRange;
        iterationImages(iterIdx) = backprojection_nearfield_mimo(iterRange, params);
        iterationMetrics(iterIdx) = imageQualityMetrics(iterationImages(iterIdx));
    end
    bestIdx = selectBestIteration(iterationMetrics);
    if bestIdx > 0 && ~isempty(iterationRanges{bestIdx})
        rangeCompensated = iterationRanges{bestIdx};
        info.bestIteration = bestIdx;
        info.bestIterationReason = 'minimum ISLR, then lower mainlobe width and entropy';
        if isfield(info.iterations(bestIdx), 'phaseTotal') && ~isempty(info.iterations(bestIdx).phaseTotal)
            phiHat = reshape(info.iterations(bestIdx).phaseTotal, [nRx 1 1 1 nPos]);
            info.phaseEstStandard = info.iterations(bestIdx).phaseTotal;
            info.phiCommonHat = squeeze(mean(info.iterations(bestIdx).phaseTotal, 1));
            info.phiRxHat = squeeze(mean(mean(info.iterations(bestIdx).phaseTotal, 2), 3));
        end
    end
end
end

function [nY, nX] = apertureGridSize(params, nPos)
nY = 1;
nX = nPos;
if isfield(params, 'sar') && isfield(params.sar, 'num_y_positions') && ...
        isfield(params.sar, 'num_x_positions') && ...
        ~isempty(params.sar.num_y_positions) && ~isempty(params.sar.num_x_positions)
    nY = params.sar.num_y_positions;
    nX = params.sar.num_x_positions;
    if nY * nX ~= nPos
        error('run_pga_mimo_compensation_experiment:BadApertureGrid', ...
            'num_y_positions * num_x_positions must match number of scan positions.');
    end
end
end

function bestIdx = selectBestIteration(iterationMetrics)
bestIdx = 0;
if isempty(iterationMetrics)
    return;
end
islr = [iterationMetrics.ISLR];
mainWidth = [iterationMetrics.mainlobeWidth];
entropy = [iterationMetrics.imageEntropy];
valid = isfinite(islr) & isfinite(mainWidth) & isfinite(entropy);
if ~any(valid)
    return;
end
score = inf(size(islr));
score(valid) = zscoreFinite(islr(valid)) + ...
    0.30 * zscoreFinite(mainWidth(valid)) + ...
    0.20 * zscoreFinite(entropy(valid));
[~, bestIdx] = min(score);
end

function z = zscoreFinite(x)
x = x(:).';
spread = std(x);
if ~isfinite(spread) || spread <= eps
    z = zeros(size(x));
else
    z = (x - mean(x)) ./ spread;
end
end

function rangeOut = applyPhaseMap(rangeIn, phiRxPos)
% phiRxPos is [rx, scanPos] in radians and is applied to complex range data.
if isempty(rangeIn.data) || any(~isfinite(rangeIn.data(:)))
    error('run_pga_mimo_compensation_experiment:BadRangeData', ...
        'rangeIn.data must be non-empty and finite.');
end
nRx = size(rangeIn.data, 2);
nPos = size(rangeIn.data, 4);
if ~isequal(size(phiRxPos), [nRx nPos])
    error('run_pga_mimo_compensation_experiment:BadPhaseMap', ...
        'phiRxPos must be [rx, scanPos].');
end
phaseFactor = reshape(exp(1j * phiRxPos), [1 nRx 1 nPos]);
rangeOut = rangeIn;
rangeOut.data = bsxfun(@times, rangeIn.data, phaseFactor);
end

function phi = makeCommonAperturePhase(nPos, amplitudeRad)
x = linspace(-1, 1, nPos);
phi = amplitudeRad * (0.80 * sin(2 * pi * 0.70 * x + 0.25) + ...
    0.30 * cos(2 * pi * 1.35 * x - 0.40));
phi = phi - mean(phi);
end

function phiRx = makeFixedRxPhase(nRx, refRx, amplitudeRad)
idx = (1:nRx).';
phiRx = amplitudeRad * sin(0.85 * idx + 0.35) + 0.25 * cos(1.70 * idx);
phiRx = phiRx - phiRx(refRx);
phiRx(refRx) = 0;
end

function [dx, dy] = makePlatformDisplacement(params, nPos)
amp = params.motion_error.amplitude_mm * 1e-3;
period = max(2, params.motion_error.period_positions);
k = (0:nPos - 1).';
phase = mod(k, period) ./ period;
triangle = 4 * abs(phase - 0.5) - 1;
dx = 0.20 * amp * sin(2 * pi * k / max(nPos - 1, 1));
dy = amp * triangle;
dx = dx - mean(dx);
dy = dy - mean(dy);
end

function phi = makeMotionEquivalentPhase(params, targets, dx, dy)
% Exact existing geometry is used for the dominant target center. This is
% equivalent to phi = 4*pi/lambda*deltaR for a monostatic first-order view,
% while the implemented MIMO form uses TX+RX path-length difference.
c = params.constants.c;
fc = params.radar.fc;
platformPos = params.sar.platform_pos;
txRel = params.array.tx_pos;
rxRel = params.array.rx_pos;
targetPos = dominantTargetPosition(targets);
nPos = size(platformPos, 1);
phi = zeros(1, nPos);

for posIdx = 1:nPos
    nominalPlatform = platformPos(posIdx, :);
    movedPlatform = nominalPlatform + [dx(posIdx), dy(posIdx), 0];
    deltaPathSum = 0;
    count = 0;
    for txIdx = 1:size(txRel, 1)
        for rxIdx = 1:size(rxRel, 1)
            txNom = nominalPlatform + txRel(txIdx, :);
            rxNom = nominalPlatform + rxRel(rxIdx, :);
            txMov = movedPlatform + txRel(txIdx, :);
            rxMov = movedPlatform + rxRel(rxIdx, :);
            nominalPath = norm(targetPos - txNom) + norm(targetPos - rxNom);
            movedPath = norm(targetPos - txMov) + norm(targetPos - rxMov);
            deltaPathSum = deltaPathSum + movedPath - nominalPath;
            count = count + 1;
        end
    end
    phi(posIdx) = -2 * pi * fc * deltaPathSum / (c * max(count, 1));
end
phi = phi - mean(phi);
end

function [nominalGradientX, nominalGradientY] = makeNominalGradients(params, targetPosition)
c = params.constants.c;
fc = params.radar.fc;
platformPos = params.sar.platform_pos;
txRel = params.array.tx_pos;
rxRel = params.array.rx_pos;
nPos = size(platformPos, 1);
[nY, nX] = apertureGridSize(params, nPos);
nTx = size(txRel, 1);
nRx = size(rxRel, 1);
signal = complex(zeros(nRx, nTx, nY, nX));

for posIdx = 1:nPos
    platform = platformPos(posIdx, :);
    [yIdx, xIdx] = ind2sub([nY nX], posIdx);
    for txIdx = 1:nTx
        txPos = platform + txRel(txIdx, :);
        for rxIdx = 1:nRx
            rxPos = platform + rxRel(rxIdx, :);
            pathLength = norm(targetPosition - txPos) + norm(targetPosition - rxPos);
            signal(rxIdx, txIdx, yIdx, xIdx) = exp(-1i * 2 * pi * fc * pathLength / c);
        end
    end
end

corrX = sum(sum(signal(:, :, :, 2:end) .* conj(signal(:, :, :, 1:end-1)), 1), 2);
corrY = sum(sum(signal(:, :, 2:end, :) .* conj(signal(:, :, 1:end-1, :)), 1), 2);
nominalGradientX = reshape(angle(corrX), [1 nY max(nX - 1, 0)]);
nominalGradientY = reshape(angle(corrY), [1 max(nY - 1, 0) nX]);
end

function weights = makeSeparableHannWeights(nY, nX)
weights = hannColumn(nY) * hannColumn(nX).';
weights = weights(:);
end

function w = hannColumn(n)
if n <= 1
    w = ones(n, 1);
else
    k = (0:n - 1).';
    w = 0.5 - 0.5 * cos(2 * pi * k / (n - 1));
end
end

function pos = dominantTargetPosition(targets)
amp = zeros(numel(targets), 1);
for idx = 1:numel(targets)
    amp(idx) = abs(targets(idx).reflectivity);
end
[~, best] = max(amp);
pos = targets(best).position;
end

function scenario = buildScenario(rangeClean, rangeDegraded, rangeCompensated, ...
    imageClean, imageDegraded, imageCompensated, phiTrue, phiCommonTrue, ...
    phiRxTrue, phiHat, info, targets)
scenario = struct();
scenario.cleanData = rangeClean.data;
scenario.degradedData = rangeDegraded.data;
scenario.compensatedData = rangeCompensated.data;
scenario.imageClean = imageClean;
scenario.imageDegraded = imageDegraded;
scenario.imageCompensated = imageCompensated;
scenario.phiTrue = phiTrue;
scenario.phiCommonTrue = phiCommonTrue;
scenario.phiRxTrue = phiRxTrue;
scenario.phiHat = phiHat;
scenario.info = info;
scenario.targets = targets;
scenario.metrics.clean = imageQualityMetrics(imageClean);
scenario.metrics.degraded = imageQualityMetrics(imageDegraded);
scenario.metrics.compensated = imageQualityMetrics(imageCompensated);
scenario.phaseCommonHat = squeeze(info.phiCommonHat);
scenario.phaseRxHat = squeeze(info.phiRxHat);
if isfield(info, 'bestIteration')
    scenario.bestIteration = info.bestIteration;
    scenario.bestIterationReason = info.bestIterationReason;
else
    scenario.bestIteration = [];
    scenario.bestIterationReason = '';
end
end

function metrics = imageQualityMetrics(image)
mag = abs(image.magnitude);
if isempty(mag) || any(~isfinite(mag(:)))
    error('run_pga_mimo_compensation_experiment:BadImage', ...
        'Image magnitude must be non-empty and finite.');
end
power = mag.^2;
[peakPower, peakIdx] = max(power(:));
[py, px] = ind2sub(size(power), peakIdx);

mainRadius = 3;
yMain = max(1, py - mainRadius):min(size(power, 1), py + mainRadius);
xMain = max(1, px - mainRadius):min(size(power, 2), px + mainRadius);
mainMask = false(size(power));
mainMask(yMain, xMain) = true;

sidePower = power(~mainMask);
if isempty(sidePower)
    sidePeak = eps;
else
    sidePeak = max(sidePower(:));
end
mainEnergy = sum(power(mainMask));
sideEnergy = sum(power(~mainMask));
prob = power(:) ./ (sum(power(:)) + eps);

metrics = struct();
metrics.PSLR = 10 * log10((sidePeak + eps) / (peakPower + eps));
metrics.ISLR = 10 * log10((sideEnergy + eps) / (mainEnergy + eps));
metrics.imageEntropy = -sum(prob .* log(prob + eps));
metrics.peakValue = sqrt(peakPower);
metrics.mainlobeWidth = estimateMainlobeWidth(power, px, py, peakPower, image);
metrics.peakPosition = image.peak_position;
metrics.energyConcentration = mainEnergy / (sum(power(:)) + eps);
end

function diagnostics = psfDiagnostics(params, imageClean)
if nargin < 2
    imageClean = [];
end
if isfield(params.scene, 'target_mode') && strcmpi(char(params.scene.target_mode), 'image')
    diagnostics = struct('name', 'all_targets', 'peaks', []);
    if isempty(imageClean)
        caseTargets = make_point_targets(params);
        caseRaw = simulate_fmcw_mimo_sar(params, caseTargets);
        caseRange = range_fft(caseRaw, params);
        imageClean = backprojection_nearfield_mimo(caseRange, params);
    end
    diagnostics.peaks = localPeakTable(imageClean, 8, 5);
    return;
end

targetTable = params.scene.targets;
numCases = size(targetTable, 1) + 1;
diagnostics = repmat(struct('name', '', 'peaks', []), numCases, 1);

for caseIdx = 1:numCases
    caseParams = params;
    if caseIdx <= size(targetTable, 1)
        caseParams.scene.targets = targetTable(caseIdx, :);
        diagnostics(caseIdx).name = sprintf('single_target_%d', caseIdx);
    else
        diagnostics(caseIdx).name = 'all_targets';
    end

    caseTargets = make_point_targets(caseParams);
    caseRaw = simulate_fmcw_mimo_sar(caseParams, caseTargets);
    caseRange = range_fft(caseRaw, caseParams);
    caseImage = backprojection_nearfield_mimo(caseRange, caseParams);
    diagnostics(caseIdx).peaks = localPeakTable(caseImage, 8, 5);
end
end

function info = targetMetadata(params, targets)
info = struct();
info.mode = 'points';
info.numScatterers = numel(targets);
if isfield(params.scene, 'target_mode')
    info.mode = char(params.scene.target_mode);
end
positions = reshape([targets.position], 3, []).';
reflectivity = [targets.reflectivity].';
info.xRange = [min(positions(:, 1)), max(positions(:, 1))];
info.yRange = [min(positions(:, 2)), max(positions(:, 2))];
info.zRange = [min(positions(:, 3)), max(positions(:, 3))];
info.maxReflectivity = max(abs(reflectivity));
if strcmpi(info.mode, 'image')
    info.imagePath = params.scene.image_path;
    info.imageNumScatterersRequested = params.scene.image_num_scatterers;
    info.imageTargetWidthM = params.scene.image_target_width_m;
    info.imageTargetCenter = params.scene.image_target_center;
end
end

function peaks = localPeakTable(image, numPeaks, suppressRadius)
mag = image.magnitude;
peakNorm = max(mag(:)) + eps;
work = mag;
peaks = zeros(numPeaks, 3);
for idx = 1:numPeaks
    [value, linearIdx] = max(work(:));
    [iy, ix] = ind2sub(size(work), linearIdx);
    peaks(idx, :) = [image.x_axis(ix), image.y_axis(iy), value / peakNorm];
    yIdx = max(1, iy - suppressRadius):min(size(work, 1), iy + suppressRadius);
    xIdx = max(1, ix - suppressRadius):min(size(work, 2), ix + suppressRadius);
    work(yIdx, xIdx) = 0;
end
end

function width = estimateMainlobeWidth(power, px, py, peakPower, image)
threshold = 0.5 * peakPower;
row = power(py, :) >= threshold;
col = power(:, px) >= threshold;
xIdx = find(row);
yIdx = find(col);
if isempty(xIdx)
    wx = 0;
else
    wx = image.x_axis(max(xIdx)) - image.x_axis(min(xIdx));
end
if isempty(yIdx)
    wy = 0;
else
    wy = image.y_axis(max(yIdx)) - image.y_axis(min(yIdx));
end
width = sqrt(wx.^2 + wy.^2);
end

function ok = isImproved(before, after)
ok = after.PSLR < before.PSLR || ...
    after.ISLR < before.ISLR || ...
    after.imageEntropy < before.imageEntropy || ...
    after.energyConcentration > before.energyConcentration;
end

function printMetricsTable(name, clean, degraded, compensated)
fprintf('\n%s\n', name);
fprintf('  Case          PSLR(dB)   ISLR(dB)   Entropy    Peak       MainWidth(m)\n');
fprintf('  clean       %9.3f %9.3f %9.4f %9.4g %12.6f\n', ...
    clean.PSLR, clean.ISLR, clean.imageEntropy, clean.peakValue, clean.mainlobeWidth);
fprintf('  degraded    %9.3f %9.3f %9.4f %9.4g %12.6f\n', ...
    degraded.PSLR, degraded.ISLR, degraded.imageEntropy, degraded.peakValue, degraded.mainlobeWidth);
fprintf('  compensated %9.3f %9.3f %9.4f %9.4g %12.6f\n', ...
    compensated.PSLR, compensated.ISLR, compensated.imageEntropy, compensated.peakValue, compensated.mainlobeWidth);
end

function printConvergenceTable(name, iterationMetrics)
if isempty(iterationMetrics)
    return;
end
fprintf('\n%s\n', name);
fprintf('  Iter          PSLR(dB)   ISLR(dB)   Entropy    Peak       MainWidth(m)  EnergyConc\n');
for idx = 1:numel(iterationMetrics)
    m = iterationMetrics(idx);
    fprintf('  %-4d       %9.3f %9.3f %9.4f %9.4g %12.6f %11.4f\n', ...
        idx, m.PSLR, m.ISLR, m.imageEntropy, m.peakValue, ...
        m.mainlobeWidth, m.energyConcentration);
end
end

function plotScenario(outDir, imageClean, scenario, figTitle, fileStem)
fig = figure('Name', figTitle, 'Color', 'w', 'Position', [100 100 1200 700]);

subplot(2, 3, 1);
plotImage(imageClean, 'Clean');
subplot(2, 3, 2);
plotImage(scenario.imageDegraded, 'Degraded');
subplot(2, 3, 3);
plotImage(scenario.imageCompensated, 'Compensated');

subplot(2, 3, 4);
plot(wrappedPhaseDeg(scenario.phiCommonTrue(:)), 'k-', 'LineWidth', 1.4, 'DisplayName', 'injected phase');
hold on;
hat = scenario.phaseCommonHat;
plot(wrappedPhaseDeg(hat(:)), 'r--', 'LineWidth', 1.2, 'DisplayName', 'estimated phase');
grid on;
xlabel('scanPos');
ylabel('phase (deg)');
title('Common aperture phase');
legend('Location', 'best');

subplot(2, 3, 5);
plot(wrappedPhaseDeg(scenario.phiRxTrue(:)), 'ko-', 'LineWidth', 1.2, 'DisplayName', 'injected phase');
hold on;
plot(wrappedPhaseDeg(scenario.phaseRxHat(:)), 'rs--', 'LineWidth', 1.2, 'DisplayName', 'estimated phase');
grid on;
xlabel('RX');
ylabel('phase (deg)');
title('Fixed RX phase');
legend('Location', 'best');

subplot(2, 3, 6);
bar([scenario.metrics.clean.PSLR, scenario.metrics.degraded.PSLR, scenario.metrics.compensated.PSLR; ...
     scenario.metrics.clean.ISLR, scenario.metrics.degraded.ISLR, scenario.metrics.compensated.ISLR; ...
     scenario.metrics.clean.imageEntropy, scenario.metrics.degraded.imageEntropy, scenario.metrics.compensated.imageEntropy].');
grid on;
set(gca, 'XTickLabel', {'clean', 'degraded', 'comp'});
legend({'PSLR', 'ISLR', 'entropy'}, 'Location', 'best');
title('Metrics');

savefig(fig, fullfile(outDir, [fileStem '.fig']));
saveas(fig, fullfile(outDir, [fileStem '.png']));
end

function plotConvergence(outDir, scenario, figTitle, fileStem)
if isempty(scenario.iterationMetrics)
    return;
end
iter = 1:numel(scenario.iterationMetrics);
pslr = [scenario.iterationMetrics.PSLR];
islr = [scenario.iterationMetrics.ISLR];
entropy = [scenario.iterationMetrics.imageEntropy];
mainWidth = [scenario.iterationMetrics.mainlobeWidth];
energyConc = [scenario.iterationMetrics.energyConcentration];
peakValue = [scenario.iterationMetrics.peakValue];

fig = figure('Name', figTitle, 'Color', 'w', 'Position', [120 120 1050 680]);

subplot(2, 3, 1);
plot(iter, pslr, 'o-', 'LineWidth', 1.3);
grid on;
xlabel('PGA iteration');
ylabel('PSLR (dB)');
title('Peak sidelobe ratio');

subplot(2, 3, 2);
plot(iter, islr, 'o-', 'LineWidth', 1.3);
grid on;
xlabel('PGA iteration');
ylabel('ISLR (dB)');
title('Integrated sidelobe ratio');

subplot(2, 3, 3);
plot(iter, entropy, 'o-', 'LineWidth', 1.3);
grid on;
xlabel('PGA iteration');
ylabel('entropy');
title('Image entropy');

subplot(2, 3, 4);
plot(iter, mainWidth, 'o-', 'LineWidth', 1.3);
grid on;
xlabel('PGA iteration');
ylabel('width (m)');
title('Mainlobe width');

subplot(2, 3, 5);
plot(iter, energyConc, 'o-', 'LineWidth', 1.3);
grid on;
xlabel('PGA iteration');
ylabel('fraction');
title('Local energy concentration');

subplot(2, 3, 6);
plot(iter, peakValue ./ max(peakValue), 'o-', 'LineWidth', 1.3);
grid on;
xlabel('PGA iteration');
ylabel('normalized peak');
title('Peak value');

savefig(fig, fullfile(outDir, [fileStem '.fig']));
saveas(fig, fullfile(outDir, [fileStem '.png']));
end

function plotImageTargetReference(outDir, targets, targetInfo)
positions = reshape([targets.position], 3, []).';
amp = abs([targets.reflectivity]).';
fig = figure('Name', 'Image target reference', 'Color', 'w', 'Position', [120 120 720 620]);
scatter(positions(:, 1), positions(:, 2), 10, amp, 'filled');
set(gca, 'YDir', 'normal');
axis image;
if exist('turbo', 'file') || exist('turbo', 'builtin')
    colormap(gca, turbo);
else
    colormap(gca, jet);
end
colorbar;
xlabel('Cross-range x (m)');
ylabel('Range y (m)');
title(sprintf('Image target reference (%d scatterers)', targetInfo.numScatterers));
savefig(fig, fullfile(outDir, 'image_target_reference.fig'));
saveas(fig, fullfile(outDir, 'image_target_reference.png'));
end

function plotImage(image, titleText)
imagesc(image.x_axis, image.y_axis, image.magnitude_db);
set(gca, 'YDir', 'normal');
axis image;
if exist('turbo', 'file') || exist('turbo', 'builtin')
    colormap(gca, turbo);
else
    colormap(gca, jet);
end
caxis([-45 0]);
colorbar;
xlabel('Cross-range x (m)');
ylabel('Range y (m)');
title(titleText);
end

function phaseDeg = wrappedPhaseDeg(phaseRad)
phaseRad = phaseRad(:);
phaseRad = phaseRad - mean(phaseRad);
phaseDeg = rad2deg(angle(exp(1i * phaseRad)));
end
