function targets = make_point_targets(params)
%MAKE_POINT_TARGETS Convert configured point targets into a structured array.

if isfield(params, 'scene') && isfield(params.scene, 'target_mode') && ...
        strcmpi(char(params.scene.target_mode), 'image')
    targets = make_image_targets(params);
    return;
end

target_table = params.scene.targets;
num_targets = size(target_table, 1);

targets(num_targets).position = [];
targets(num_targets).reflectivity = [];

for idx = 1:num_targets
    targets(idx).position = real(target_table(idx, 1:3));
    targets(idx).reflectivity = target_table(idx, 4);
end
end

function targets = make_image_targets(params)
image_path = params.scene.image_path;
resolved_path = resolve_image_path(image_path);
if isempty(resolved_path)
    error('make_point_targets:ImageNotFound', ...
        'params.scene.image_path must point to an existing image file.');
end

[img, ~, alpha] = imread(resolved_path);
img = im2double(img);
if ndims(img) == 3
    gray = 0.2989 * img(:, :, 1) + 0.5870 * img(:, :, 2) + 0.1140 * img(:, :, 3);
else
    gray = img;
end

amplitude = 1 - gray;
if ~isempty(alpha)
    amplitude = amplitude .* im2double(alpha);
end
amplitude = amplitude - min(amplitude(:));
if max(amplitude(:)) > eps
    amplitude = amplitude ./ max(amplitude(:));
end

num_scatterers = get_scene_scalar(params, 'image_num_scatterers', 3000);
valid = amplitude > 0;
if ~any(valid(:))
    error('make_point_targets:EmptyImageTarget', ...
        'Image target contains no non-background scatterers.');
end

[amp_sorted, ~] = sort(amplitude(:), 'descend');
amp_sorted = amp_sorted(amp_sorted > 0);
n_keep = min(num_scatterers, numel(amp_sorted));
cutoff = amp_sorted(n_keep);
order = find(amplitude(:) >= cutoff);
if numel(order) > n_keep
    sample_idx = unique(round(linspace(1, numel(order), n_keep)));
    while numel(sample_idx) < n_keep
        missing = setdiff(1:numel(order), sample_idx, 'stable');
        sample_idx = [sample_idx(:); missing(1:n_keep - numel(sample_idx)).']; %#ok<AGROW>
        sample_idx = sort(sample_idx);
    end
    order = order(sample_idx(1:n_keep));
end
[row_idx, col_idx] = ind2sub(size(amplitude), order);
amp = amplitude(order);

target_width = get_scene_scalar(params, 'image_target_width_m', 0.060);
target_center = get_scene_vector(params, 'image_target_center', default_image_center(params), 2);
target_height = target_width * size(amplitude, 1) / size(amplitude, 2);

if size(amplitude, 2) == 1
    x = target_center(1) * ones(size(col_idx));
else
    x = target_center(1) + ((col_idx - 1) / (size(amplitude, 2) - 1) - 0.5) * target_width;
end
if size(amplitude, 1) == 1
    y = target_center(2) * ones(size(row_idx));
else
    y = target_center(2) + (0.5 - (row_idx - 1) / (size(amplitude, 1) - 1)) * target_height;
end
z = params.scene.z0 * ones(size(x));

targets(n_keep).position = [];
targets(n_keep).reflectivity = [];
for idx = 1:n_keep
    targets(idx).position = [x(idx), y(idx), z(idx)];
    targets(idx).reflectivity = amp(idx);
end
end

function resolved_path = resolve_image_path(image_path)
resolved_path = '';
if ~(ischar(image_path) || isstring(image_path))
    return;
end

candidate = char(image_path);
if exist(candidate, 'file')
    resolved_path = candidate;
    return;
end

if numel(candidate) >= 7 && strcmp(candidate(1:5), '/mnt/') && candidate(7) == '/'
    drive = upper(candidate(6));
    rest = strrep(candidate(8:end), '/', filesep);
    windows_candidate = [drive, ':', filesep, rest];
    if exist(windows_candidate, 'file')
        resolved_path = windows_candidate;
    end
end
end

function value = get_scene_scalar(params, name, default_value)
value = default_value;
if isfield(params.scene, name) && ~isempty(params.scene.(name))
    value = params.scene.(name);
end
end

function value = get_scene_vector(params, name, default_value, n)
value = default_value;
if isfield(params.scene, name) && ~isempty(params.scene.(name))
    value = params.scene.(name);
end
value = reshape(value, 1, []);
if numel(value) ~= n
    error('make_point_targets:BadSceneVector', ...
        'params.scene.%s must have %d elements.', name, n);
end
end

function center = default_image_center(params)
center = [mean(params.image.x_axis), mean(params.image.y_axis)];
end
