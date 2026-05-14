function [dataCorr, phaseEst, info] = phaseCompensateMIMO(data, opts)
%PHASECOMPENSATEMIMO Estimate and compensate position-dependent RX phase jitter.
%
%   [dataCorr, phaseEst, info] = phaseCompensateMIMO(data, opts)
%
%   Default data layout is data(rx, tx, freq, y, x). The default method is a
%   PGA-like phase-gradient estimator on the 2-D scan aperture. The older
%   reference-channel coherence estimator is still available through
%   opts.method = 'coherence'.
%
%   opts is optional. Supported fields:
%       refRx                  Reference RX index. Default: 1
%       dimOrder               'rx_tx_freq_y_x' or a 5-token permutation.
%                              Default: 'rx_tx_freq_y_x'
%       freqIdx                Frequency indices used for estimation.
%                              Default: all frequencies
%       method                 'pga', 'rx_phase', 'pga_plus_rx', or
%                              'coherence'. Default: 'pga'
%       numIterations          Number of PGA compensation iterations.
%                              Default: 1
%       pgaRangeBins           Range bins used by PGA after IFFT over freq.
%                              Default: strongest bins are selected
%       pgaNumRangeBins        Number of strongest range bins if pgaRangeBins
%                              is empty. Default: 8
%       pgaInputDomain         'frequency' or 'range'. Use 'range' when the
%                              third dimension is already range-compressed.
%                              Default: 'frequency'
%       phaseModel             'rx_relative', 'aperture_common', or
%                              'mimo_aperture'. Default: 'rx_relative'
%       nominalGradientX       Optional nominal aperture phase gradient in x.
%       nominalGradientY       Optional nominal aperture phase gradient in y.
%       smoothWindow           [y x] moving-average window for phase maps.
%                              Default: [5 5]
%       unwrapSpatial          Spatially unwrap phase maps. Default: true
%       removeMeanPhase        Remove each RX phase-map weighted mean.
%                              Default: true
%       removeLinearPhase      Remove best-fit x/y linear phase ramps after
%                              estimation. Default: false
%       minMagnitudePercentile Low-strength positions are downweighted.
%                              Default: 20
%       enableAmplitudeEq      Normalize each RX by its mean magnitude
%                              before phase estimation. Default: false
%       restoreAmplitudeAfterEq Restore RX amplitudes before returning
%                              corrected data. Default: true
%
%   phaseEst is returned in the original data dimension order with singleton
%   dimensions at tx and freq axes, so it can be inspected or broadcast.
if nargin < 2 || isempty(opts)
    opts = struct();
end
validateattributes(data, {'numeric'}, {}, mfilename, 'data', 1);
opts = mergeOptions(opts);

if ~isfloat(data)
    data = double(data);
end
if isempty(data) || any(~isfinite(data(:)))
    error('phaseCompensateMIMO:BadData', ...
        'data must be non-empty and contain only finite values.');
end

[dataStd, standardPerm, dimNames] = toStandardOrder(data, opts.dimOrder);
sz = size(dataStd);
if numel(sz) < 5
    sz(end+1:5) = 1;
end
nRx = sz(1);
nFreq = sz(3);
nY = sz(4);
nX = sz(5);

if opts.refRx > nRx
    error('phaseCompensateMIMO:BadReferenceRx', ...
        'opts.refRx=%d exceeds number of RX channels (%d).', opts.refRx, nRx);
end

freqIdx = opts.freqIdx;
if isempty(freqIdx)
    freqIdx = 1:nFreq;
end
if any(freqIdx > nFreq)
    error('phaseCompensateMIMO:BadFrequencyIndex', ...
        'opts.freqIdx contains an index greater than number of frequencies (%d).', nFreq);
end

[dataStd, gainRx] = amplitudeEqualizeRx(dataStd, opts.enableAmplitudeEq);

switch lower(opts.method)
    case 'pga'
        [dataCorrStd, phaseStd, weightStd, maskStd, methodInfo, iterationInfo] = ...
            compensatePGAIterative(dataStd, freqIdx, opts);
    case 'rx_phase'
        [phaseStd, weightStd, maskStd, methodInfo] = estimateFixedRxPhase(dataStd, freqIdx, opts);
        comp = exp(-1i * reshape(phaseStd, [nRx 1 1 nY nX]));
        dataCorrStd = dataStd .* comp;
        iterationInfo = [];
    case 'pga_plus_rx'
        commonOpts = opts;
        commonOpts.method = 'pga';
        commonOpts.phaseModel = 'aperture_common';
        [dataCommonStd, phaseCommonStd, weightStd, maskStd, methodInfoCommon, iterationInfo] = ...
            compensatePGAIterative(dataStd, freqIdx, commonOpts);
        [phaseRxStd, weightRxStd, maskRxStd, methodInfoRx] = ...
            estimateFixedRxPhase(dataCommonStd, freqIdx, opts);
        compRx = exp(-1i * reshape(phaseRxStd, [nRx 1 1 nY nX]));
        dataCorrStd = dataCommonStd .* compRx;
        phaseStd = phaseCommonStd + phaseRxStd;
        weightStd = max(weightStd, weightRxStd);
        maskStd = maskStd & maskRxStd;
        methodInfo = methodInfoCommon;
        methodInfo.description = 'Common aperture PGA followed by fixed RX phase calibration';
        methodInfo.phaseCommon = squeeze(phaseCommonStd(1, :, :));
        methodInfo.phaseRxFixed = squeeze(phaseRxStd(:, 1, 1));
        methodInfo.rxPhase = methodInfoRx;
    case 'coherence'
        [phaseStd, weightStd, maskStd, methodInfo] = estimatePhaseCoherence(dataStd, freqIdx, opts);
        comp = exp(-1i * reshape(phaseStd, [nRx 1 1 nY nX]));
        dataCorrStd = dataStd .* comp;
        iterationInfo = [];
    otherwise
        error('phaseCompensateMIMO:BadMethod', ...
            'opts.method must be ''pga'', ''rx_phase'', ''pga_plus_rx'', or ''coherence''.');
end

phaseBroadcastStd = reshape(phaseStd, [nRx 1 1 nY nX]);
if opts.enableAmplitudeEq && opts.restoreAmplitudeAfterEq
    dataCorrStd = restoreRxAmplitude(dataCorrStd, gainRx);
end
dataCorr = ipermute(dataCorrStd, standardPerm);
phaseEst = ipermute(phaseBroadcastStd, standardPerm);

info = struct();
info.dimOrder = opts.dimOrder;
info.standardDimOrder = 'rx_tx_freq_y_x';
info.originalDimNames = dimNames;
info.method = lower(opts.method);
info.refRx = opts.refRx;
info.freqIdx = freqIdx;
info.phaseEstStandard = phaseStd;
info.correlationMagnitude = weightStd;
info.validMask = maskStd;
info.methodInfo = methodInfo;
info.iterations = iterationInfo;
info.coherenceBefore = channelCoherence(dataStd, opts.refRx, freqIdx);
info.coherenceAfter = channelCoherence(dataCorrStd, opts.refRx, freqIdx);
info.gainRx = gainRx;
info.phiCommonHat = extractCommonPhase(methodInfo, phaseStd);
info.phiRxHat = extractRxPhase(methodInfo, phaseStd, opts.refRx);
info.phiRxHatByAperture = phaseStd;
end

function opts = mergeOptions(opts)
if ~isstruct(opts)
    error('phaseCompensateMIMO:BadOptions', 'opts must be a structure.');
end

if isfield(opts, 'numIter') && ~isfield(opts, 'numIterations')
    opts.numIterations = opts.numIter;
end
if isfield(opts, 'selectedRangeBins') && ~isfield(opts, 'pgaRangeBins')
    opts.pgaRangeBins = opts.selectedRangeBins;
end
if isfield(opts, 'unwrapPhase') && ~isfield(opts, 'unwrapSpatial')
    opts.unwrapSpatial = opts.unwrapPhase;
end

defaults = struct();
defaults.refRx = 1;
defaults.dimOrder = 'rx_tx_freq_y_x';
defaults.freqIdx = [];
defaults.method = 'pga';
defaults.numIterations = 1;
defaults.pgaRangeBins = [];
defaults.pgaNumRangeBins = 8;
defaults.pgaInputDomain = 'frequency';
defaults.phaseModel = 'rx_relative';
defaults.nominalGradientX = [];
defaults.nominalGradientY = [];
defaults.smoothWindow = [5 5];
defaults.unwrapSpatial = true;
defaults.removeMeanPhase = true;
defaults.removeLinearPhase = false;
defaults.minMagnitudePercentile = 20;
defaults.enableAmplitudeEq = false;
defaults.restoreAmplitudeAfterEq = true;
defaults.selectedImageRegion = [];

names = fieldnames(defaults);
for k = 1:numel(names)
    name = names{k};
    if ~isfield(opts, name) || isempty(opts.(name))
        opts.(name) = defaults.(name);
    end
end

validateattributes(opts.refRx, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, mfilename, 'opts.refRx');
validateattributes(opts.numIterations, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, mfilename, 'opts.numIterations');
if ~isempty(opts.freqIdx)
    validateattributes(opts.freqIdx, {'numeric'}, ...
        {'vector', 'integer', 'positive'}, mfilename, 'opts.freqIdx');
end
if ~isempty(opts.pgaRangeBins)
    validateattributes(opts.pgaRangeBins, {'numeric'}, ...
        {'vector', 'integer', 'positive'}, mfilename, 'opts.pgaRangeBins');
end
validateattributes(opts.pgaNumRangeBins, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, mfilename, 'opts.pgaNumRangeBins');
validateattributes(opts.smoothWindow, {'numeric'}, ...
    {'numel', 2, 'integer', 'positive'}, mfilename, 'opts.smoothWindow');
validateattributes(opts.minMagnitudePercentile, {'numeric'}, ...
    {'scalar', '>=', 0, '<=', 100}, mfilename, 'opts.minMagnitudePercentile');
if ~(ischar(opts.method) || isstring(opts.method))
    error('phaseCompensateMIMO:BadMethod', 'opts.method must be text.');
end
if ~(ischar(opts.pgaInputDomain) || isstring(opts.pgaInputDomain))
    error('phaseCompensateMIMO:BadPgaInputDomain', ...
        'opts.pgaInputDomain must be ''frequency'' or ''range''.');
end
if ~(ischar(opts.phaseModel) || isstring(opts.phaseModel))
    error('phaseCompensateMIMO:BadPhaseModel', ...
        'opts.phaseModel must be ''rx_relative'', ''aperture_common'', or ''mimo_aperture''.');
end

opts.method = char(opts.method);
opts.numIterations = double(opts.numIterations);
opts.pgaInputDomain = lower(char(opts.pgaInputDomain));
if ~ismember(opts.pgaInputDomain, {'frequency', 'range'})
    error('phaseCompensateMIMO:BadPgaInputDomain', ...
        'opts.pgaInputDomain must be ''frequency'' or ''range''.');
end
opts.phaseModel = lower(char(opts.phaseModel));
if ~ismember(opts.phaseModel, {'rx_relative', 'aperture_common', 'mimo_aperture'})
    error('phaseCompensateMIMO:BadPhaseModel', ...
        'opts.phaseModel must be ''rx_relative'', ''aperture_common'', or ''mimo_aperture''.');
end
opts.smoothWindow = reshape(double(opts.smoothWindow), 1, 2);
opts.unwrapSpatial = logical(opts.unwrapSpatial);
opts.removeMeanPhase = logical(opts.removeMeanPhase);
opts.removeLinearPhase = logical(opts.removeLinearPhase);
opts.enableAmplitudeEq = logical(opts.enableAmplitudeEq);
opts.restoreAmplitudeAfterEq = logical(opts.restoreAmplitudeAfterEq);
end

function [dataStd, gainRx] = amplitudeEqualizeRx(dataStd, enabled)
nRx = size(dataStd, 1);
gainRx = ones(nRx, 1);
if ~enabled
    return;
end

for rx = 1:nRx
    rxData = dataStd(rx, :, :, :, :);
    gain = mean(abs(rxData(:)));
    if ~isfinite(gain) || gain <= eps
        gain = 1;
    end
    gainRx(rx) = gain;
    dataStd(rx, :, :, :, :) = dataStd(rx, :, :, :, :) ./ gain;
end
end

function dataStd = restoreRxAmplitude(dataStd, gainRx)
for rx = 1:size(dataStd, 1)
    dataStd(rx, :, :, :, :) = dataStd(rx, :, :, :, :) .* gainRx(rx);
end
end

function [phaseStd, weightStd, maskStd, methodInfo] = estimateFixedRxPhase(dataStd, freqIdx, opts)
nRx = size(dataStd, 1);
nY = size(dataStd, 4);
nX = size(dataStd, 5);
estData = dataStd(:, :, freqIdx, :, :);
refData = estData(opts.refRx, :, :, :, :);

phaseRx = zeros(nRx, 1);
weightRx = zeros(nRx, 1);
for rx = 1:nRx
    cur = estData(rx, :, :, :, :);
    corr = sum(cur(:) .* conj(refData(:)));
    phaseRx(rx) = angle(corr);
    weightRx(rx) = abs(corr);
end
phaseRx(opts.refRx) = 0;

phaseStd = repmat(reshape(phaseRx, [nRx 1 1]), [1 nY nX]);
weightStd = repmat(reshape(weightRx, [nRx 1 1]), [1 nY nX]);
maskStd = buildMagnitudeMask(weightStd, opts.minMagnitudePercentile);

methodInfo = struct();
methodInfo.description = 'Fixed RX phase calibration relative to reference RX';
methodInfo.gradientX = [];
methodInfo.gradientY = [];
methodInfo.rangeBins = [];
methodInfo.phaseCommon = zeros(1, nY, nX);
methodInfo.phaseRxRelative = phaseStd;
methodInfo.phaseRxFixed = phaseRx;
end

function phiCommon = extractCommonPhase(methodInfo, phaseStd)
if isfield(methodInfo, 'phaseCommon') && ~isempty(methodInfo.phaseCommon)
    phiCommon = squeeze(methodInfo.phaseCommon);
elseif size(phaseStd, 1) >= 1
    phiCommon = squeeze(mean(phaseStd, 1));
else
    phiCommon = [];
end
end

function phiRx = extractRxPhase(methodInfo, phaseStd, refRx)
if isfield(methodInfo, 'phaseRxFixed') && ~isempty(methodInfo.phaseRxFixed)
    phiRx = methodInfo.phaseRxFixed(:);
elseif isfield(methodInfo, 'phaseRxRelative') && ~isempty(methodInfo.phaseRxRelative)
    tmp = methodInfo.phaseRxRelative;
    phiRx = squeeze(mean(mean(tmp, 2), 3));
else
    phiRx = squeeze(mean(mean(phaseStd, 2), 3));
end
if ~isempty(phiRx) && refRx <= numel(phiRx)
    phiRx = angle(exp(1i * (phiRx - phiRx(refRx))));
end
end

function [phaseStd, weightStd, maskStd, methodInfo] = estimatePhaseCoherence(dataStd, freqIdx, opts)
nRx = size(dataStd, 1);
nY = size(dataStd, 4);
nX = size(dataStd, 5);

estData = dataStd(:, :, freqIdx, :, :);
refData = estData(opts.refRx, :, :, :, :);
refData = repmat(refData, [nRx 1 1 1 1]);

corrMap = reshape(sum(sum(estData .* conj(refData), 2), 3), [nRx nY nX]);

phaseStd = angle(corrMap);
weightStd = abs(corrMap);
phaseStd(opts.refRx, :, :) = 0;

maskStd = buildMagnitudeMask(weightStd, opts.minMagnitudePercentile);
phaseStd = smoothPhaseMaps(phaseStd, weightStd, maskStd, opts);
phaseStd(opts.refRx, :, :) = 0;

methodInfo = struct();
methodInfo.description = 'Reference RX coherence phase estimate';
methodInfo.gradientX = [];
methodInfo.gradientY = [];
methodInfo.rangeBins = [];
end

function [dataCorrStd, phaseTotal, weightStd, maskStd, methodInfo, iterationInfo] = ...
    compensatePGAIterative(dataStd, freqIdx, opts)
nRx = size(dataStd, 1);
nY = size(dataStd, 4);
nX = size(dataStd, 5);

dataIter = dataStd;
phaseTotal = zeros(nRx, nY, nX);
iterationInfo = repmat(struct( ...
    'phaseStep', [], ...
    'phaseTotal', [], ...
    'dataStd', [], ...
    'coherenceBefore', [], ...
    'coherenceAfter', [], ...
    'methodInfo', []), opts.numIterations, 1);

for iterIdx = 1:opts.numIterations
    [phaseStep, weightStd, maskStd, methodInfo] = estimatePhasePGA(dataIter, freqIdx, opts);
    comp = exp(-1i * reshape(phaseStep, [nRx 1 1 nY nX]));
    coherenceBefore = channelCoherence(dataIter, opts.refRx, freqIdx);
    dataIter = dataIter .* comp;
    phaseTotal = phaseTotal + phaseStep;

    iterationInfo(iterIdx).phaseStep = phaseStep;
    iterationInfo(iterIdx).phaseTotal = phaseTotal;
    iterationInfo(iterIdx).dataStd = dataIter;
    iterationInfo(iterIdx).coherenceBefore = coherenceBefore;
    iterationInfo(iterIdx).coherenceAfter = channelCoherence(dataIter, opts.refRx, freqIdx);
    iterationInfo(iterIdx).methodInfo = methodInfo;
end

dataCorrStd = dataIter;
methodInfo.numIterations = opts.numIterations;
end

function [phaseStd, weightStd, maskStd, methodInfo] = estimatePhasePGA(dataStd, freqIdx, opts)
nRx = size(dataStd, 1);
nFreq = numel(freqIdx);
nY = size(dataStd, 4);
nX = size(dataStd, 5);

estData = dataStd(:, :, freqIdx, :, :);
switch opts.pgaInputDomain
    case 'frequency'
        rangeData = ifft(estData, [], 3);
    case 'range'
        rangeData = estData;
    otherwise
        error('phaseCompensateMIMO:BadPgaInputDomain', ...
            'opts.pgaInputDomain must be ''frequency'' or ''range''.');
end
rangeBins = selectPgaRangeBins(rangeData, opts, nFreq);
pgaData = rangeData(:, :, rangeBins, :, :);

prodX = pgaData(:, :, :, :, 2:end) .* conj(pgaData(:, :, :, :, 1:end-1));
prodY = pgaData(:, :, :, 2:end, :) .* conj(pgaData(:, :, :, 1:end-1, :));

corrX = reshape(sum(sum(prodX, 2), 3), [nRx nY max(nX - 1, 0)]);
corrY = reshape(sum(sum(prodY, 2), 3), [nRx max(nY - 1, 0) nX]);

gradXRaw = angle(corrX);
gradYRaw = angle(corrY);
weightX = abs(corrX);
weightY = abs(corrY);

[nominalGradX, nominalGradY] = nominalGradients(opts, nY, nX);
[phaseStd, phaseCommon, phaseRxRelative, gradX, gradY] = ...
    estimatePhaseByModel(corrX, corrY, gradXRaw, gradYRaw, weightX, weightY, ...
    nominalGradX, nominalGradY, nRx, nY, nX, opts);

weightStd = apertureWeightFromGradients(weightX, weightY, nRx, nY, nX);
maskStd = buildMagnitudeMask(weightStd, opts.minMagnitudePercentile);
phaseStd = smoothPhaseMaps(phaseStd, weightStd, maskStd, opts);
if opts.removeLinearPhase
    phaseStd = removeLinearPhaseMaps(phaseStd, weightStd, maskStd);
end
if strcmp(opts.phaseModel, 'rx_relative')
    phaseStd(opts.refRx, :, :) = 0;
end

methodInfo = struct();
methodInfo.description = 'PGA-like 2-D aperture phase-gradient estimate';
methodInfo.pgaInputDomain = opts.pgaInputDomain;
methodInfo.phaseModel = opts.phaseModel;
methodInfo.rangeBins = rangeBins;
methodInfo.gradientX = gradX;
methodInfo.gradientY = gradY;
methodInfo.gradientWeightX = weightX;
methodInfo.gradientWeightY = weightY;
methodInfo.nominalGradientX = nominalGradX;
methodInfo.nominalGradientY = nominalGradY;
methodInfo.phaseCommon = phaseCommon;
methodInfo.phaseRxRelative = phaseRxRelative;
end

function [phaseStd, phaseCommon, phaseRxRelative, gradX, gradY] = ...
    estimatePhaseByModel(corrX, corrY, gradXRaw, gradYRaw, weightX, weightY, ...
    nominalGradX, nominalGradY, nRx, nY, nX, opts)
phaseCommon = zeros(1, nY, nX);
phaseRxRelative = zeros(nRx, nY, nX);
gradX = gradXRaw;
gradY = gradYRaw;

if any(strcmp(opts.phaseModel, {'aperture_common', 'mimo_aperture'}))
    corrCommonX = reshape(sum(corrX, 1), [1 nY max(nX - 1, 0)]);
    corrCommonY = reshape(sum(corrY, 1), [1 max(nY - 1, 0) nX]);
    gradCommonX = angle(corrCommonX);
    gradCommonY = angle(corrCommonY);

    if nX > 1
        gradCommonX = angle(exp(1i * (gradCommonX - nominalGradX)));
    end
    if nY > 1
        gradCommonY = angle(exp(1i * (gradCommonY - nominalGradY)));
    end

    weightCommonX = reshape(sum(abs(corrX), 1), [1 nY max(nX - 1, 0)]);
    weightCommonY = reshape(sum(abs(corrY), 1), [1 max(nY - 1, 0) nX]);
    phaseCommon = integratePgaGradients( ...
        gradCommonX, gradCommonY, weightCommonX, weightCommonY, 1, nY, nX);
    phaseCommon = removeWeightedMeanPhase(phaseCommon, ...
        apertureWeightFromGradients(weightCommonX, weightCommonY, 1, nY, nX));
    gradX = repmat(gradCommonX, [nRx 1 1]);
    gradY = repmat(gradCommonY, [nRx 1 1]);
end

if any(strcmp(opts.phaseModel, {'rx_relative', 'mimo_aperture'}))
    gradRxX = gradXRaw;
    gradRxY = gradYRaw;
    refGradX = gradRxX(opts.refRx, :, :);
    refGradY = gradRxY(opts.refRx, :, :);
    gradRxX = angle(exp(1i * (gradRxX - repmat(refGradX, [nRx 1 1]))));
    gradRxY = angle(exp(1i * (gradRxY - repmat(refGradY, [nRx 1 1]))));
    phaseRxRelative = integratePgaGradients(gradRxX, gradRxY, weightX, weightY, nRx, nY, nX);
    phaseRxRelative(opts.refRx, :, :) = 0;
    gradX = gradRxX;
    gradY = gradRxY;
end

switch opts.phaseModel
    case 'rx_relative'
        phaseStd = phaseRxRelative;
    case 'aperture_common'
        phaseStd = repmat(phaseCommon, [nRx 1 1]);
    case 'mimo_aperture'
        phaseStd = repmat(phaseCommon, [nRx 1 1]) + phaseRxRelative;
    otherwise
        error('phaseCompensateMIMO:BadPhaseModel', 'Unsupported phase model.');
end
end

function [nominalGradX, nominalGradY] = nominalGradients(opts, nY, nX)
nominalGradX = zeros(1, nY, max(nX - 1, 0));
nominalGradY = zeros(1, max(nY - 1, 0), nX);

if ~isempty(opts.nominalGradientX)
    nominalGradX = reshape(opts.nominalGradientX, [1 nY max(nX - 1, 0)]);
end
if ~isempty(opts.nominalGradientY)
    nominalGradY = reshape(opts.nominalGradientY, [1 max(nY - 1, 0) nX]);
end
end

function phaseOut = removeWeightedMeanPhase(phaseIn, weight)
phaseOut = phaseIn;
for idx = 1:size(phaseIn, 1)
    phaseMap = reshape(phaseIn(idx, :, :), size(phaseIn, 2), size(phaseIn, 3));
    weightMap = reshape(weight(idx, :, :), size(weight, 2), size(weight, 3));
    denom = sum(weightMap(:));
    if denom > 0
        phaseMap = phaseMap - sum(phaseMap(:) .* weightMap(:)) / denom;
    end
    phaseOut(idx, :, :) = phaseMap;
end
end

function rangeBins = selectPgaRangeBins(rangeData, opts, nFreq)
if ~isempty(opts.pgaRangeBins)
    if any(opts.pgaRangeBins > nFreq)
        error('phaseCompensateMIMO:BadPgaRangeBins', ...
            'opts.pgaRangeBins contains an index greater than selected frequency/range length (%d).', nFreq);
    end
    rangeBins = opts.pgaRangeBins(:).';
    return;
end

rangePower = sum(abs(rangeData).^2, 1);
rangePower = sum(rangePower, 2);
rangePower = sum(rangePower, 4);
rangePower = sum(rangePower, 5);
rangePower = reshape(rangePower, [1 nFreq]);

nKeep = min(opts.pgaNumRangeBins, nFreq);
[~, order] = sort(rangePower, 'descend');
rangeBins = sort(order(1:nKeep));
end

function phaseStd = integratePgaGradients(gradX, gradY, weightX, weightY, nRx, nY, nX)
phaseStd = zeros(nRx, nY, nX);

for rx = 1:nRx
    phaseMap = zeros(nY, nX);
    support = false(nY, nX);
    support(1, 1) = true;

    for x = 2:nX
        phaseMap(1, x) = phaseMap(1, x - 1) + gradX(rx, 1, x - 1);
        support(1, x) = true;
    end
    for y = 2:nY
        phaseMap(y, 1) = phaseMap(y - 1, 1) + gradY(rx, y - 1, 1);
        support(y, 1) = true;
    end

    for y = 2:nY
        for x = 2:nX
            wx = weightX(rx, y, x - 1);
            wy = weightY(rx, y - 1, x);
            fromX = phaseMap(y, x - 1) + gradX(rx, y, x - 1);
            fromY = phaseMap(y - 1, x) + gradY(rx, y - 1, x);
            denom = wx + wy;

            if denom > eps
                phaseMap(y, x) = (wx * fromX + wy * fromY) / denom;
            else
                phaseMap(y, x) = 0.5 * (fromX + fromY);
            end
            support(y, x) = true;
        end
    end

    phaseMap(~support) = 0;
    phaseStd(rx, :, :) = phaseMap;
end
end

function weightStd = apertureWeightFromGradients(weightX, weightY, nRx, nY, nX)
weightStd = zeros(nRx, nY, nX);

if nX > 1
    weightStd(:, :, 1:end-1) = weightStd(:, :, 1:end-1) + weightX;
    weightStd(:, :, 2:end) = weightStd(:, :, 2:end) + weightX;
end

if nY > 1
    weightStd(:, 1:end-1, :) = weightStd(:, 1:end-1, :) + weightY;
    weightStd(:, 2:end, :) = weightStd(:, 2:end, :) + weightY;
end

for rx = 1:nRx
    plane = reshape(weightStd(rx, :, :), nY, nX);
    if all(plane(:) == 0)
        plane(:) = 1;
    end
    weightStd(rx, :, :) = plane;
end
end

function [dataStd, standardPerm, dimNames] = toStandardOrder(data, dimOrder)
if ~(ischar(dimOrder) || isstring(dimOrder))
    error('phaseCompensateMIMO:BadDimOrder', 'opts.dimOrder must be text.');
end

dimNames = strsplit(char(dimOrder), '_');
if numel(dimNames) ~= 5
    error('phaseCompensateMIMO:BadDimOrder', ...
        'opts.dimOrder must contain five tokens: rx, tx, freq, y, x.');
end

required = {'rx', 'tx', 'freq', 'y', 'x'};
if ~all(ismember(required, dimNames)) || numel(unique(dimNames)) ~= 5
    error('phaseCompensateMIMO:BadDimOrder', ...
        'opts.dimOrder must be a permutation of rx_tx_freq_y_x.');
end

standardPerm = zeros(1, 5);
for k = 1:5
    standardPerm(k) = find(strcmp(dimNames, required{k}), 1);
end

dataStd = permute(data, standardPerm);
end

function mask = buildMagnitudeMask(weight, percentileValue)
mask = true(size(weight));
if percentileValue <= 0
    return;
end

for rx = 1:size(weight, 1)
    plane = reshape(weight(rx, :, :), size(weight, 2), size(weight, 3));
    threshold = prctile(plane(:), percentileValue);
    mask(rx, :, :) = plane >= threshold;
end
end

function phaseOut = smoothPhaseMaps(phaseIn, weight, mask, opts)
phaseOut = phaseIn;
win = opts.smoothWindow;
kernel = ones(win(1), win(2));

for rx = 1:size(phaseIn, 1)
    phaseMap = reshape(phaseIn(rx, :, :), size(phaseIn, 2), size(phaseIn, 3));
    weightMap = reshape(weight(rx, :, :), size(weight, 2), size(weight, 3));
    maskMap = reshape(mask(rx, :, :), size(mask, 2), size(mask, 3));

    if opts.unwrapSpatial
        phaseMap = unwrap(unwrap(phaseMap, [], 1), [], 2);
    end

    weightMap(~maskMap) = 0;
    if opts.removeMeanPhase
        denom = sum(weightMap(:));
        if denom > 0
            phaseMap = phaseMap - sum(phaseMap(:) .* weightMap(:)) / denom;
        end
    end

    if any(win > 1)
        numerator = conv2(phaseMap .* weightMap, kernel, 'same');
        denominator = conv2(weightMap, kernel, 'same');
        smoothMap = phaseMap;
        valid = denominator > eps;
        smoothMap(valid) = numerator(valid) ./ denominator(valid);
        phaseMap = smoothMap;
    end

    phaseOut(rx, :, :) = phaseMap;
end
end

function phaseOut = removeLinearPhaseMaps(phaseIn, weight, mask)
phaseOut = phaseIn;
nY = size(phaseIn, 2);
nX = size(phaseIn, 3);
[xGrid, yGrid] = meshgrid(1:nX, 1:nY);

for rx = 1:size(phaseIn, 1)
    phaseMap = reshape(phaseIn(rx, :, :), nY, nX);
    weightMap = reshape(weight(rx, :, :), nY, nX);
    maskMap = reshape(mask(rx, :, :), nY, nX);
    valid = maskMap & isfinite(phaseMap) & weightMap > 0;

    if nnz(valid) < 3
        continue;
    end

    xValid = xGrid(valid);
    yValid = yGrid(valid);
    design = ones(nnz(valid), 1);
    if nX > 1
        design = [design, xValid(:)];
    end
    if nY > 1
        design = [design, yValid(:)];
    end
    response = phaseMap(valid);
    response = response(:);
    w = sqrt(weightMap(valid));
    w = w(:);
    coeff = bsxfun(@times, design, w) \ (response .* w);
    ramp = coeff(1) * ones(nY, nX);
    coeffIdx = 2;
    if nX > 1
        ramp = ramp + coeff(coeffIdx) * xGrid;
        coeffIdx = coeffIdx + 1;
    end
    if nY > 1
        ramp = ramp + coeff(coeffIdx) * yGrid;
    end
    phaseOut(rx, :, :) = phaseMap - ramp;
end
end

function coh = channelCoherence(dataStd, refRx, freqIdx)
nRx = size(dataStd, 1);
coh = zeros(nRx, 1);
ref = dataStd(refRx, :, freqIdx, :, :);

for rx = 1:nRx
    cur = dataStd(rx, :, freqIdx, :, :);
    numerator = abs(sum(cur(:) .* conj(ref(:))));
    denominator = sqrt(sum(abs(cur(:)).^2) * sum(abs(ref(:)).^2));
    if denominator > 0
        coh(rx) = numerator / denominator;
    else
        coh(rx) = NaN;
    end
end
end
