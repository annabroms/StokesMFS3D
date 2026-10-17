function grid = build_image_discretization(q,opt)
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
% (empty selects the paper's gap rule), imageFamilies (subset of SRTD),
% capMode, capAngles, capUpsampling or exact capPointCounts, capLayers,
% capExponent, capMinAzimuth,
% rowWeightMode, and useImages. grid.body(k) owns its basis blocks,
% collocation nodes, physical area weights, image lines, and pair metadata.
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
assert(opt.capMode ~= "exponential" || opt.rowWeightMode == "local", ...
    'Exponentially clustered caps require local merged-grid row areas.');
P = size(q,1);
assert(P >= 1,'At least one sphere is required.');
[proxy,coarseCollocation] = getDesignGrid(opt.proxyRadius, ...
    struct('des_n',opt.N,'a_glob',opt.collocationRatio));
template = struct('family',"",'points',zeros(0,3), ...
    'normalDirections',zeros(0,3),'scale',1);
bodyTemplate = struct('basis',repmat(template,0,1), ...
    'collocation',zeros(0,3),'weights',zeros(0,1), ...
    'imagePoints',zeros(0,3),'imageDirections',zeros(0,3), ...
    'capPoints',zeros(0,3),'contactDirections',zeros(0,3), ...
    'N',size(proxy,1),'N_image',0);
bodies = repmat(bodyTemplate,P,1);
for k = 1:P
    bodies(k).collocation = q(k,:)+coarseCollocation;
    bodies(k).weights = (4*pi/size(coarseCollocation,1))* ...
        ones(size(coarseCollocation,1),1);
end
pairs = zeros(0,4); % [first body, second body, gap, image nodes]
for k = 1:P-1
    for j = k+1:P
        gap = norm(q(j,:)-q(k,:))-2;
        assert(gap > 0,'Spheres must not overlap or touch.');
        if gap >= opt.nearThreshold
            continue
        end
        if isempty(opt.N_image)
            N_image = max(0,ceil(-8.72*log10(gap)-6.15));
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
            if N_image == 0
                % Keep gap metadata for independent validation, but leave
                % the basic source and collocation grids untouched.
                continue
            end
            if opt.useImages
                [line,~] = image_line_points(q(body,:),q(neighbor,:), ...
                    1,opt.proxyRadius,N_image,opt.imageOffsetFraction);
                bodies(body).imagePoints = [bodies(body).imagePoints;line];
                bodies(body).imageDirections = [bodies(body).imageDirections; ...
                    repmat(direction,N_image,1)];
                bodies(body).N_image = bodies(body).N_image+N_image;
            end
            sphere = struct('axes',ones(3,1),'Q',eye(3),'c',q(body,:).');
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
                    sphere,q(body,:).'+direction.',capOptions);
                capArea = 2*pi*(1-cos(opt.capAngles(cap)))/capCount;
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
        sphere = struct('axes',ones(3,1),'Q',eye(3),'c',q(k,:).');
        bodies(k).weights = collocation_patch_area( ...
            bodies(k).collocation,ones(size(bodies(k).collocation,1),1), ...
            {sphere});
    end
    block = template;
    block.family = "S";
    block.points = q(k,:)+proxy;
    bodies(k).basis = block;
    for family = char(opt.imageFamilies)
        if bodies(k).N_image == 0 || ~opt.useImages
            break
        end
        block = template;
        block.family = string(family);
        block.points = bodies(k).imagePoints;
        if family == 'T'
            block.normalDirections = bodies(k).imageDirections;
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
    'rowWeightMode',"local",'useImages',true,'imageOffsetFraction',0.05);
q = [0,0,0;2.01,0,0;4.02,0,0];
grid = build_image_discretization(q,opt);
assert(isequal([grid.body.N_image],[3,6,3]));
assert(size(grid.body(2).collocation,1) > ...
    size(grid.body(1).collocation,1));
for k = 1:3
    assert(all(grid.body(k).weights > 0));
    assert(abs(sum(grid.body(k).weights)-4*pi) < 1e-9);
end
opt.imageFamilies = "SRTD";
opt.capPointCounts = [24,16];
stressletGrid = build_image_discretization(q,opt);
assert(numel(stressletGrid.body(1).basis) == 5);
assert(size(stressletGrid.body(1).capPoints,1) == 40);
assert(isequal(size(stressletGrid.body(1).basis(4).normalDirections), ...
    size(stressletGrid.body(1).imageDirections)));
fprintf('build_image_discretization self-test passed.\n');
end
