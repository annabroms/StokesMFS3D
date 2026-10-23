function grid = build_image_discretization(q,opt,radii)
%BUILD_IMAGE_DISCRETIZATION Assemble sphere-owned image MFS nodes and areas.
%   grid = build_image_discretization(q,opt)
%   build_image_discretization() runs a small unequal-contact self-test.
% Reference: A. Broms, A. H. Barnett, and A.-K. Tornberg,
% "Accurate close interactions of Stokes spheres using lubrication-adapted
% image systems," Journal of Computational Physics 523 (2025), 113636.
%
% q is P-by-3 of nonoverlapping unit-sphere centres. opt is the completed
% option set from solve_resistance_with_images; important controls are
% proxyRadius, N, collocationRatio, nearThreshold, N_image
% (empty selects the paper's gap rule), imageFamilies
% (S/R/D/T/T_t/T_n/T_n2/Q_s),
% capMode, capAngles, capUpsampling or exact capPointCounts, capLayers,
% capExponent, capMinAzimuth,
% rowWeightMode, useRowScaling, and useImages. grid.body(k) owns its basis
% blocks, collocation nodes, physical area weights, effective row scales,
% image lines, and pair metadata.
% Every vector family has interleaved [x1;y1;z1;x2;...] coefficients.
% The coarse design carries Stokeslets regardless of imageFamilies.
%
% The reference cap weights are 4*pi/M0 on coarse nodes and
% 2*pi*(1-cos(beta))/Mcap on each Fibonacci cap. "local" weights instead
% use positive lumped areas on the deduplicated complete surface grid.

if nargin == 0
    self_test();
    grid = struct();
    return
end
validateattributes(q,{'numeric'},{'real','finite','2d','ncols',3});
if nargin < 3 || isempty(radii)
    radii = ones(size(q,1),1);
end
radii = radii(:);
validateattributes(radii,{'numeric'}, ...
    {'real','finite','positive','numel',size(q,1)});
if ~isfield(opt,'useRowScaling')
    opt.useRowScaling = true; % Also support direct callers with older opt structs.
end
if ~isfield(opt,'use_orth_dir_images')
    opt.use_orth_dir_images = false;
end
assert(opt.capMode ~= "exponential" || opt.rowWeightMode == "local", ...
    'Exponentially clustered caps require local merged-grid row areas.');
P = size(q,1);
families = image_family_list(opt.imageFamilies);
assert(P >= 1,'At least one sphere is required.');
[proxy,coarseCollocation] = getDesignGrid(opt.proxyRadius, ...
    struct('des_n',opt.N,'a_glob',opt.collocationRatio));
template = struct('family',"",'points',zeros(0,3), ...
    'normalDirections',zeros(0,3),'scale',1);
bodyTemplate = struct('basis',repmat(template,0,1), ...
    'collocation',zeros(0,3),'weights',zeros(0,1), ...
    'rowScale',zeros(0,1), ...
    'imagePoints',zeros(0,3),'imageDirections',zeros(0,3), ...
    'capPoints',zeros(0,3),'contactDirections',zeros(0,3), ...
    'N',size(proxy,1),'N_image',0);
bodies = repmat(bodyTemplate,P,1);
for k = 1:P
    bodies(k).collocation = q(k,:)+radii(k)*coarseCollocation;
    bodies(k).weights = (4*pi*radii(k)^2/size(coarseCollocation,1))* ...
        ones(size(coarseCollocation,1),1);
end
pairs = zeros(0,4); % [first body, second body, gap, image nodes]
for k = 1:P-1
    for j = k+1:P
        gap = norm(q(j,:)-q(k,:))-radii(k)-radii(j);
        assert(gap > 0,'Spheres must not overlap or touch.');
        relativeGap = gap/min(radii(k),radii(j));
        if relativeGap >= opt.nearThreshold
            continue
        end
        if isempty(opt.N_image)
            N_image = max(0,ceil(-8.72*log10(relativeGap)-6.15));
        else
            N_image = opt.N_image;
        end
        pairs(end+1,:) = [k,j,gap,N_image]; %#ok<AGROW>
        for side = 1:2
            if side == 1
                body = k;
                neighbor = j;
            else
                body = j;
                neighbor = k;
            end
            direction = (q(neighbor,:)-q(body,:))/norm(q(neighbor,:)-q(body,:));
            bodies(body).contactDirections(end+1,:) = direction;
            if N_image == 0 || ~opt.useImages
                % Keep gap metadata for independent validation, but leave
                % the basic source and collocation grids untouched.
                continue
            end
            [line,~] = image_line_points(q(body,:),q(neighbor,:), ...
                radii(body),radii(body)*opt.proxyRadius,N_image, ...
                opt.imageOffsetFraction,radii(neighbor));
            bodies(body).imagePoints = [bodies(body).imagePoints;line];
            bodies(body).imageDirections = [bodies(body).imageDirections; ...
                repmat(direction,N_image,1)];
            bodies(body).N_image = bodies(body).N_image+N_image;
            sphere = struct('axes',radii(body)*ones(3,1), ...
                'Q',eye(3),'c',q(body,:).');
            % Exponential grading already refines toward contact, so one
            % cap holds the points of both fibonacci caps.
            nCaps = 2;
            if opt.capMode == "exponential"
                nCaps = 1;
            end
            for cap = 1:nCaps
                capCount = round(opt.capUpsampling(cap)*3*N_image);
                if isfield(opt,'capPointCounts') && ...
                        ~isempty(opt.capPointCounts)
                    capCount = opt.capPointCounts(cap);
                end
                if nCaps == 1
                    capCount = round(sum(opt.capUpsampling)*3*N_image);
                    if isfield(opt,'capPointCounts') && ...
                            ~isempty(opt.capPointCounts)
                        capCount = sum(opt.capPointCounts);
                    end
                end
                capOptions.collocationPatchMode = opt.capMode;
                capOptions.collocationPatchAngle = opt.capAngles(cap);
                capOptions.collocationPatchPointCount = capCount;
                minAz = 3;
                if opt.capMode == "exponential"
                    minAz = opt.capMinAzimuth;
                end
                % Honour the minimum points per ring when the image rule
                % gives fewer cap nodes than the requested layer count needs.
                capOptions.collocationPatchLayers = min( ...
                    opt.capLayers(cap),floor((capCount-1)/minAz));
                assert(capOptions.collocationPatchLayers >= 1, ...
                    'Each cap needs at least four collocation nodes.');
                capOptions.collocationExponent = opt.capExponent(cap);
                capOptions.collocationMinAzimuth = opt.capMinAzimuth;
                capOptions.collocationAzimuthOversampling = 1;
                [capPoints,~] = collocation_surface_patch( ...
                    sphere,q(body,:).'+radii(body)*direction.',capOptions);
                capArea = 2*pi*radii(body)^2* ...
                    (1-cos(opt.capAngles(cap)))/capCount;
                for node = 1:size(capPoints,1)
                    point = capPoints(node,:);
                    % Keep the first occurrence, including coarse design nodes.
                    if any(sum((bodies(body).collocation-point).^2,2) < 1e-20)
                        continue
                    end
                    bodies(body).collocation(end+1,:) = point;
                    bodies(body).capPoints(end+1,:) = point;
                    bodies(body).weights(end+1,1) = capArea;
                end
            end
        end
    end
end

for k = 1:P
    if opt.rowWeightMode == "local"
        sphere = struct('axes',radii(k)*ones(3,1), ...
            'Q',eye(3),'c',q(k,:).');
        bodies(k).weights = collocation_patch_area( ...
            bodies(k).collocation,ones(size(bodies(k).collocation,1),1), ...
            {sphere});
    end
    if opt.useRowScaling
        bodies(k).rowScale = sqrt(bodies(k).weights);
    else
        bodies(k).rowScale = ones(size(bodies(k).weights));
    end
    block = template;
    block.family = "S";
    block.points = q(k,:)+radii(k)*proxy;
    bodies(k).basis = block;
    for family = families
        if bodies(k).N_image == 0 || ~opt.useImages
            break
        end
        block = template;
        block.family = string(family);
        block.points = bodies(k).imagePoints;
        if ismember(family,["T","T_n","T_n2","T_t"])
            block.normalDirections = image_stresslet_directions( ...
                bodies(k).imageDirections,opt.use_orth_dir_images,family);
        elseif family=="Q_s"
            block.normalDirections = bodies(k).imageDirections;
            block.scale = radii(k)^2;
        end
        bodies(k).basis(end+1,1) = block;
    end
    numberOfColumns = 3*sum(arrayfun( ...
        @(b) size(b.points,1),bodies(k).basis));
    assert(3*size(bodies(k).collocation,1) > numberOfColumns, ...
        'Each one-body collocation system must be overdetermined.');
end
grid.body = bodies;
grid.centres = q;
grid.radii = radii;
grid.pairs = pairs;
grid.N = size(proxy,1);
grid.coarseCollocationCount = size(coarseCollocation,1);
end

function self_test()
%SELF_TEST checks that contact counts and body-owned blocks can differ.
root = fileparts(fileparts(mfilename('fullpath')));
oldPath = path; cleanup = onCleanup(@() path(oldPath));
addpath(fullfile(root,'geometry','spheres'), ...
    fullfile(root,'geometry','ellipsoids'));
opt = struct('proxyRadius',0.63,'N',50,'collocationRatio',2, ...
    'nearThreshold',0.15,'N_image',3,'imageFamilies',"SRD", ...
    'capMode',"fibonacci",'capAngles',[pi/5,pi/60], ...
    'capUpsampling',[6,6],'capLayers',[4,4], ...
    'capExponent',[3,3],'capMinAzimuth',3, ...
    'rowWeightMode',"local",'useRowScaling',true, ...
    'useImages',true,'imageOffsetFraction',0.05);
q = [0,0,0;2.01,0,0;4.02,0,0];
grid = build_image_discretization(q,opt);
assert(isequal([grid.body.N_image],[3,6,3]));
assert(size(grid.body(2).collocation,1) > ...
    size(grid.body(1).collocation,1));
for k = 1:3
    assert(all(grid.body(k).weights > 0));
    assert(abs(sum(grid.body(k).weights)-4*pi) < 1e-9);
    assert(norm(grid.body(k).rowScale-sqrt(grid.body(k).weights)) == 0);
end
opt.useRowScaling = false;
plainGrid = build_image_discretization(q,opt);
for k = 1:3
    assert(isequal(plainGrid.body(k).collocation,grid.body(k).collocation));
    assert(isequal(plainGrid.body(k).weights,grid.body(k).weights));
    assert(all(plainGrid.body(k).rowScale == 1));
end
opt.useRowScaling = true;
opt.imageFamilies = "SRTD";
opt.capPointCounts = [24,16];
stressletGrid = build_image_discretization(q,opt);
assert(numel(stressletGrid.body(1).basis) == 5);
assert(size(stressletGrid.body(1).capPoints,1) == 40);
assert(isequal(size(stressletGrid.body(1).basis(4).normalDirections), ...
    size(stressletGrid.body(1).imageDirections)));
assert(isequal(stressletGrid.body(1).basis(4).normalDirections, ...
    stressletGrid.body(1).imageDirections));
opt.use_orth_dir_images = true;
orthogonalGrid = build_image_discretization(q,opt);
for k = 1:3
    ray = orthogonalGrid.body(k).imageDirections;
    direction = orthogonalGrid.body(k).basis(4).normalDirections;
    assert(max(abs(sum(ray.*direction,2)))<1e-14);
    assert(isequal(orthogonalGrid.body(k).imagePoints, ...
        stressletGrid.body(k).imagePoints));
end
opt.useImages = false;
basicGrid = build_image_discretization(q,opt);
assert(all(arrayfun(@(body) isscalar(body.basis), basicGrid.body)));
assert(all(arrayfun(@(body) isempty(body.capPoints), basicGrid.body)));
assert(all(arrayfun(@(body) isempty(body.imagePoints), basicGrid.body)));
assert(~isempty(basicGrid.pairs));
fprintf('build_image_discretization self-test passed.\n');
end
