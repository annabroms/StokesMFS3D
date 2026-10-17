function [basicGrid,enhancedGrid] = lap3d_build_two_body_grids(pair,opt,basicGrid)
%LAP3D_BUILD_TWO_BODY_GRIDS Build basic and locally enhanced MFS grids.
%
% pair must contain E1, E2, p1, p2, n, Pminus, Pplus, and gap.  Each body
% structure contains axes, rotation Q, and centre c.  The two bodies must
% have the same semi-axes because the basic grid uses one reference body.
% Request only one output to build the basic grid alone.  A sweep may pass
% a cached basicGrid as the third input when only the enhancing sources change.
% Equal spheres use the continuous proxy radius and mirrored enhancing caps;
% the existing sampled-proxy enclosure rule is preserved for ellipsoids.
% opt.removeDuplicates controls source deduplication only.  When false,
% retain all chopped source nodes; collocation deduplication is unchanged.
% opt.sourceChopOffset adds a signed physical depth to the proxy cut:
% positive extends the retained cap inside the basic proxy, negative stops
% short.  This does not change the parent geometry or its sampling.
% Basic nodes always carry Laplace monopoles.  SLP+DLP adds one normal
% dipole on each retained cap node; normals point out of its parent surface.
% Dipole indices refer to the unchanged body-ordered physical source array.
% Optional Laplace preconditioning metadata preserves physical grid locations.
% Local surface areas are estimated on complete merged grids, only if enabled.

requiredFields = {'E1','E2','p1','p2','n','Pminus','Pplus','gap'};
for k = 1:numel(requiredFields)
    if ~isfield(pair,requiredFields{k})
        error('lap3d_build_two_body_grids:MissingPairField', ...
            'pair.%s is required.',requiredFields{k});
    end
end

bodies = {pair.E1,pair.E2};
if norm(pair.E1.axes-pair.E2.axes,inf) > ...
        100*eps*max(1,norm(pair.E1.axes,inf))
    error('lap3d_build_two_body_grids:DifferentBodyAxes', ...
        'The repository basic-grid generator requires matching body axes.');
end

if nargin < 3 || isempty(basicGrid)
    basicGrid = build_basic_grid(pair,bodies,opt);
end
if nargout > 1
    enhancedGrid = build_enhanced_grid(pair,bodies,basicGrid,opt);
end
end


function grid = build_basic_grid(pair,bodies,opt)
%BUILD_BASIC_GRID maps repository nodes to inward proxies and body surfaces.

semiAxes = pair.E1.axes(:).';
rotations = {pair.E1.Q,pair.E2.Q};
centres = [pair.E1.c.';pair.E2.c.'];

if opt.basicGridDistribution == "sphericaldesign"
    designOptions.des_n = opt.basicSourceCount;
    designOptions.a_glob = opt.basicCollocationCount/opt.basicSourceCount;
    [sourceDesign,collocationDesign] = getDesignGrid(1,designOptions);

    % Preserve inward NORMAL displacement, also for mapped ellipsoid nodes.
    sourceSurface = sourceDesign.*semiAxes;
    sourceNormals = sourceDesign./semiAxes;
    sourceNormals = sourceNormals./vecnorm(sourceNormals,2,2);
    referenceSources = sourceSurface-opt.sourceOffset*sourceNormals;
    referenceCollocation = collocationDesign.*semiAxes;
    sources = zeros(2*size(referenceSources,1),3);
    collocation = zeros(2*size(referenceCollocation,1),3);
    for body = 1:2
        sourceRows = (body-1)*size(referenceSources,1) + ...
            (1:size(referenceSources,1));
        collocationRows = (body-1)*size(referenceCollocation,1) + ...
            (1:size(referenceCollocation,1));
        sources(sourceRows,:) = ...
            (bodies{body}.c+bodies{body}.Q*referenceSources.').';
        collocation(collocationRows,:) = ...
            (bodies{body}.c+bodies{body}.Q*referenceCollocation.').';
    end
else
    [sources,collocation] = getEllipsoidGrids( ...
        semiAxes,2,pair.gap,opt.coarseSourceResolution(1), ...
        opt.coarseSourceResolution(2),opt.sourceOffset,rotations,centres, ...
        opt.basicCollocationResolution);
end

sourceCount = size(sources,1)/2;
collocationCount = size(collocation,1)/2;
if sourceCount ~= round(sourceCount) || collocationCount ~= round(collocationCount)
    error('lap3d_build_two_body_grids:UnequalBodyBlocks', ...
        'The repository grid did not produce equal-sized body blocks.');
end

sourceBlocks = {sources(1:sourceCount,:),sources(sourceCount+1:end,:)};
collocationBlocks = {collocation(1:collocationCount,:), ...
    collocation(collocationCount+1:end,:)};
[testPoints,testBody,testNearContact,testMetadata] = lap3d_evaluation_grid(pair,opt);

grid.sources = sources;
grid.basicGridDistribution = opt.basicGridDistribution;
grid.collocation = collocation;
grid.sourceBody = [ones(sourceCount,1);2*ones(sourceCount,1)];
grid.collocationBody = [ones(collocationCount,1); ...
    2*ones(collocationCount,1)];
grid.testPoints = testPoints;
grid.testBody = testBody;
grid.testNearContact = testNearContact;
grid.testCounts = testMetadata.counts;
grid.testPatchCounts = testMetadata.patchCounts;
grid.testGridMetadata = testMetadata;
grid.sourceBlocks = sourceBlocks;
grid.collocationBlocks = collocationBlocks;
grid.sourceCapBlocks = {zeros(0,3),zeros(0,3)};
grid.enhancingSourceMask = false(size(sources,1),1);
grid.enhancingLaplaceBasis = "SLP";
grid.dipoleSourceIndices = zeros(0,1);
grid.dipoleNormals = zeros(0,3);
grid.dipoleLengthScale = opt.laplaceDipoleLengthScale;
grid.columnScaling = opt.columnScaling;
grid.surfaceElementPreconditioning = opt.surfaceElementPreconditioning;
grid.collocationSurfaceElements = [];
if opt.surfaceElementPreconditioning
    grid.collocationSurfaceElements = collocation_patch_area( ...
        grid.collocation,grid.collocationBody,bodies);
end
grid.collocationPatchBlocks = {zeros(0,3),zeros(0,3)};
grid.enhancingCollocationCounts = [0;0];
grid.totalEnhancingCollocation = 0;
grid.focalDiagnostic = focal_segment_diagnostic(pair.E1,opt.sourceOffset);

validate_grid_points(grid,bodies);
end


function diagnostic = focal_segment_diagnostic(E,sourceOffset)
%FOCAL_SEGMENT_DIAGNOSTIC checks the prolate or spherical focal set.

semiAxes = E.axes(:);
axisTolerance = 1e-12*max(semiAxes);
% if abs(semiAxes(2)-semiAxes(3)) > axisTolerance || ...
%         semiAxes(1) < semiAxes(2)-axisTolerance
%     error('lap3d_build_two_body_grids:NotProlate', ...
%         'The focal diagnostic requires axes [a,b,b] with a>=b.');
% end

if abs(semiAxes(1)-semiAxes(2)) <= axisTolerance
    diagnostic.focalHalfLength = 0;
    diagnostic.proxyIsConfocal = true;
else
    diagnostic.focalHalfLength = ...
        sqrt(semiAxes(1)^2-semiAxes(2)^2);
    diagnostic.proxyIsConfocal = false;
end

diagnostic.focalEndpoints = [ ...
    E.c-E.Q(:,1)*diagnostic.focalHalfLength, ...
    E.c+E.Q(:,1)*diagnostic.focalHalfLength];
diagnostic.tipClearance = semiAxes(1)-diagnostic.focalHalfLength;
diagnostic.enclosureMargin = diagnostic.tipClearance-sourceOffset;
diagnostic.enclosed = diagnostic.enclosureMargin > 0;

if ~diagnostic.enclosed
    warning('lap3d_build_two_body_grids:FocalSegmentOutside', ...
        ['The inward normal-offset source surface does not enclose the ', ...
         'complete focal segment.']);
end
end


function grid = build_enhanced_grid(pair,bodies,basicGrid,opt)
%BUILD_ENHANCED_GRID appends chopped source caps and graded target patches.

contacts = {pair.p1,pair.p2};
accumulation = {pair.Pminus,pair.Pplus};
outwardDirections = {pair.n,-pair.n};

sourceBlocks = cell(2,1);
collocationBlocks = cell(2,1);
sourceCapBlocks = cell(2,1);
collocationPatchBlocks = cell(2,1);
sourceCapMetadata = cell(2,1);
collocationPatchMetadata = cell(2,1);
equalSpheres = max(pair.E1.axes)-min(pair.E1.axes) < ...
    1e-12*max(pair.E1.axes);
reflectionCentre = 0.5*(pair.E1.c+pair.E2.c);
reflectionNormal = pair.E2.c-pair.E1.c;
reflectionNormal = reflectionNormal/norm(reflectionNormal);
reflection = eye(3)-2*(reflectionNormal*reflectionNormal.');

for body = 1:2
    if equalSpheres && body == 2
        % Reflect the retained cap, including any duplicates if requested.
        sourceCap = (reflectionCentre+ ...
            reflection*(sourceCapBlocks{1}.'-reflectionCentre)).';
        metadata = sourceCapMetadata{1};
        metadata.parentPoints = (reflectionCentre+ ...
            reflection*(metadata.parentPoints.'-reflectionCentre)).';
        pointFields = {'parentCentre','nearFocus','innerFocus', ...
            'nearAccumulationPoint','radialSurfacePoint','outerTip', ...
            'nearestBasicProxyPoint','proxyReferencePoint'};
        for k = 1:numel(pointFields)
            name = pointFields{k};
            metadata.(name) = reflectionCentre+ ...
                reflection*(metadata.(name)-reflectionCentre);
        end
        % Flip a transverse axis too so the stored frame remains proper.
        metadata.parentRotation = reflection*metadata.parentRotation* ...
            diag([1,-1,1]);
        metadata.isMirrored = true;
        sourceCapMetadata{body} = metadata;
    else
        [sourceCap,sourceCapMetadata{body}] = build_source_cap( ...
            bodies{body},contacts{body},accumulation{body}, ...
            outwardDirections{body},basicGrid.sourceBlocks{body},opt);
        if opt.removeDuplicates
            sourceCap = remove_near_duplicates(sourceCap, ...
                basicGrid.sourceBlocks{body},1e-11*max(bodies{body}.axes));
        end
    end

    [collocationPatch,collocationPatchMetadata{body}] = ...
        collocation_surface_patch(bodies{body},contacts{body},opt);
    collocationPatch = remove_near_duplicates(collocationPatch, ...
        basicGrid.collocationBlocks{body},1e-11*max(bodies{body}.axes));
    collocationPatchMetadata{body}.retainedCount = size(collocationPatch,1);

    sourceCapMetadata{body}.retainedCount = size(sourceCap,1);
    sourceCapBlocks{body} = sourceCap;
    collocationPatchBlocks{body} = collocationPatch;
    sourceBlocks{body} = [basicGrid.sourceBlocks{body};sourceCap];
    collocationBlocks{body} = [basicGrid.collocationBlocks{body}; ...
        collocationPatch];
end

grid = basicGrid;
grid.sources = vertcat(sourceBlocks{:});
grid.collocation = vertcat(collocationBlocks{:});
grid.sourceBody = [ones(size(sourceBlocks{1},1),1); ...
    2*ones(size(sourceBlocks{2},1),1)];
grid.collocationBody = [ones(size(collocationBlocks{1},1),1); ...
    2*ones(size(collocationBlocks{2},1),1)];
grid.sourceBlocks = sourceBlocks;
grid.collocationBlocks = collocationBlocks;
grid.sourceCapBlocks = sourceCapBlocks;
grid.collocationPatchBlocks = collocationPatchBlocks;
grid.sourceCapMetadata = sourceCapMetadata;
grid.collocationPatchMetadata = collocationPatchMetadata;
grid.enhancingCollocationCounts = [size(collocationPatchBlocks{1},1); ...
    size(collocationPatchBlocks{2},1)];
grid.totalEnhancingCollocation = sum(grid.enhancingCollocationCounts);
grid.columnScaling = opt.columnScaling;
grid.surfaceElementPreconditioning = opt.surfaceElementPreconditioning;
grid.collocationSurfaceElements = [];
if opt.surfaceElementPreconditioning
    % Recompute on the MERGED grid; do not add separate overlapping grid areas.
    grid.collocationSurfaceElements = collocation_patch_area( ...
        grid.collocation,grid.collocationBody,bodies);
end

grid.enhancingSourceMask = [false(size(basicGrid.sourceBlocks{1},1),1); ...
    true(size(sourceCapBlocks{1},1),1); ...
    false(size(basicGrid.sourceBlocks{2},1),1); ...
    true(size(sourceCapBlocks{2},1),1)];
grid.enhancingLaplaceBasis = opt.enhancingLaplaceBasis;
grid.dipoleSourceIndices = zeros(0,1);
grid.dipoleNormals = zeros(0,3);
grid.dipoleLengthScale = opt.laplaceDipoleLengthScale;
if opt.enhancingLaplaceBasis == "SLP+DLP"
    normalBlocks = cell(2,1);
    for body = 1:2
        parent = sourceCapMetadata{body};
        parentQuadratic = parent.parentRotation* ...
            diag(1./parent.parentAxes.^2)*parent.parentRotation.';
        normalBlocks{body} = ...
            (sourceCapBlocks{body}-parent.parentCentre.')*parentQuadratic;
        normalBlocks{body} = normalBlocks{body} ./ ...
            vecnorm(normalBlocks{body},2,2);
    end
    grid.dipoleSourceIndices = find(grid.enhancingSourceMask);
    grid.dipoleNormals = vertcat(normalBlocks{:});
end

if size(grid.collocation,1) <= size(grid.sources,1)
    error('lap3d_build_two_body_grids:TooFewCollocationNodes', ...
        ['The enhanced grid has %d collocation nodes and %d sources. ', ...
         'Increase collocationPatchLayers or collocationMinAzimuth, or ', ...
         'decrease enhancingSourceCount (design) or ', ...
         'fineSourceResolution (quadrature).'], ...
        size(grid.collocation,1),size(grid.sources,1));
end

validate_grid_points(grid,bodies);
end


function [cap,metadata] = build_source_cap(E,contactPoint, ...
        accumulationPoints,outwardDirection,basicProxyPoints,opt)
%BUILD_SOURCE_CAP makes a chopped spheroid or sphere around the focal set.

outwardDirection = outwardDirection/norm(outwardDirection);
isSphere = max(E.axes)-min(E.axes) < 1e-12*max(E.axes);
proxyRadius = NaN;
proxyReferenceIsSampled = ~isSphere;
accumulationDistances = vecnorm(accumulationPoints-contactPoint,2,1);
[~,nearestIndex] = min(accumulationDistances);
nearAccumulationPoint = accumulationPoints(:,nearestIndex);

accumulationPhysicalLevel = ellipsoid_level(E,accumulationPoints.');
if max(accumulationPhysicalLevel) >= 1
    error('lap3d_build_two_body_grids:AccumulationOutsideBody', ...
        'Every accumulation point must lie strictly inside its body.');
end

if opt.enhancingSourceShape == "spheroid"
    % As in ellipse_cheb_segment, the body centre and nearest accumulation
    % point are the foci before the contact-facing tip is stretched.
    focalAxis = nearAccumulationPoint-E.c;
    focalDistance = norm(focalAxis);
    if focalDistance == 0
        error('lap3d_build_two_body_grids:DegenerateFocalAxis', ...
            'The accumulation point must differ from the body centre.');
    end
    focalAxis = focalAxis/focalDistance;

    localAxis = E.Q.'*focalAxis;
    bodyRadialDistance = 1/sqrt(sum((localAxis./E.axes).^2));
    radialSurfacePoint = E.c+bodyRadialDistance*focalAxis;
    distanceToSurface = bodyRadialDistance-focalDistance;
    if distanceToSurface <= 0
        error('lap3d_build_two_body_grids:BadAccumulationPoint', ...
            'The nearest accumulation point must lie inside the radial surface.');
    end

    innerFocus = E.c;
    nearFocus = nearAccumulationPoint;
    focalHalfLength = focalDistance/2;
    proxyCentre = 0.5*(innerFocus+nearFocus);
    proxyRotation = frame_with_z_axis(focalAxis,E.Q(:,1));

    requestedTipExtension = opt.sourceBeta*distanceToSurface;
    majorRadius = focalHalfLength+requestedTipExtension;

    if isSphere
        % Use the continuous proxy radius, independent of sampled sources.
        proxyRadius = E.axes(1)-opt.sourceOffset;
        proxyReferencePoint = E.c+proxyRadius*focalAxis;
        nearestBasicProxyPoint = NaN(3,1);  % no sampled-point enclosure test
    else
        proxyDistances = vecnorm(basicProxyPoints-nearAccumulationPoint.',2,2);
        [~,nearestProxyIndex] = min(proxyDistances);
        nearestBasicProxyPoint = basicProxyPoints(nearestProxyIndex,:).';
        proxyReferencePoint = nearestBasicProxyPoint;
    end
    localProxyPoint = proxyRotation.'*(proxyReferencePoint-proxyCentre);
    axialCoordinate = localProxyPoint(3);
    % Sphere bounds are axial and independent of width.  Ellipsoid bounds
    % include the transverse scaling when enclosing the sampled proxy point.
    if isSphere
        requiredMajorRadius = abs(proxyRadius-focalHalfLength)*(1+1e-10);
    else
        transverseRadiusSquared = ...
            sum(localProxyPoint(1:2).^2)/opt.sourceWidthScale^2;
        quadraticCoefficient = focalHalfLength^2 + ...
            axialCoordinate^2 + transverseRadiusSquared;
        discriminant = quadraticCoefficient^2 - ...
            4*axialCoordinate^2*focalHalfLength^2;
        requiredMajorSquared = 0.5*(quadraticCoefficient + ...
            sqrt(max(0,discriminant)));
        requiredMajorRadius = sqrt(requiredMajorSquared)*(1+1e-10);
    end
    majorRadius = max(majorRadius,requiredMajorRadius);
    tipExtension = majorRadius-focalHalfLength;
    minorRadius = opt.sourceWidthScale* ...
        sqrt(majorRadius^2-focalHalfLength^2);
    if minorRadius > majorRadius
        error('lap3d_build_two_body_grids:OblateSourceParent', ...
            'sourceWidthScale makes the parent oblate; reduce the width scale.');
    end
else
    innerFocus = nearAccumulationPoint;
    nearFocus = nearAccumulationPoint;
    focalHalfLength = 0;
    tipExtension = opt.sourceSphereRadius;
    majorRadius = opt.sourceSphereRadius;
    minorRadius = opt.sourceSphereRadius;
    proxyCentre = nearAccumulationPoint;
    proxyRotation = frame_with_z_axis(outwardDirection,E.Q(:,1));
    radialSurfacePoint = contactPoint;
    distanceToSurface = min(accumulationDistances);
    requestedTipExtension = tipExtension;
    if isSphere
        proxyRadius = E.axes(1)-opt.sourceOffset;
        proxyReferencePoint = E.c+proxyRadius*outwardDirection;
        nearestBasicProxyPoint = NaN(3,1);
    else
        proxyDistances = vecnorm(basicProxyPoints-nearAccumulationPoint.',2,2);
        [~,nearestProxyIndex] = min(proxyDistances);
        nearestBasicProxyPoint = basicProxyPoints(nearestProxyIndex,:).';
        proxyReferencePoint = nearestBasicProxyPoint;
    end
    requiredMajorRadius = NaN;
end

if opt.enhancingSourceDistribution == "sphericaldesign"
    % Map a complete unit-sphere design onto the unchanged parent surface.
    sourceDesign = get_sphdesign(opt.enhancingSourceCount);
    parentNodes = sourceDesign.*[minorRadius,minorRadius,majorRadius];
else
    parent = ellipsoid_param(minorRadius,minorRadius,majorRadius);
    sourceQuadratureOptions.tensor = ...
        opt.enhancingSourceDistribution == "tensor";
    sourceQuadratureOptions.extraunodes = ...
        opt.enhancingSourceExtraAzimuthNodes;
    sourceQuadratureOptions.minunodes = opt.enhancingSourceMinAzimuthNodes;
    parent = setupsurfquad(parent,opt.fineSourceResolution, ...
        sourceQuadratureOptions);
    parentNodes = parent.x.';
end
parentPoints = (proxyCentre+proxyRotation*parentNodes.').';
physicalLevel = ellipsoid_level(E,parentPoints);

% For an aligned prolate parent in a sphere, the outer tip has the greatest
% distance to the body centre.  This checks the whole continuous parent.
if isSphere && norm(proxyCentre-E.c)+majorRadius >= E.axes(1)
    error('lap3d_build_two_body_grids:SourceParentOutsideBody', ...
        'The parent outer tip leaves the sphere; reduce sourceBeta or sourceSphereRadius.');
end
if max(physicalLevel) >= 1
    error('lap3d_build_two_body_grids:SourceParentOutsideBody', ...
        ['The parent source surface leaves the physical body (max level %.6g). ', ...
         'Reduce sourceBeta, sourceWidthScale, or sourceSphereRadius.'], ...
        max(physicalLevel));
end

proxyCoordinates = proxyRotation.'*(accumulationPoints-proxyCentre);
accumulationLevel = (proxyCoordinates(1,:)/minorRadius).^2 + ...
    (proxyCoordinates(2,:)/minorRadius).^2 + ...
    (proxyCoordinates(3,:)/majorRadius).^2;
if max(accumulationLevel) >= 1
    error('lap3d_build_two_body_grids:AccumulationNotEnclosed', ...
        'The parent source surface does not enclose every accumulation point.');
end

nearestProxyCoordinates = proxyRotation.'*(proxyReferencePoint-proxyCentre);
nearestProxyLevel = sum((nearestProxyCoordinates./ ...
    [minorRadius;minorRadius;majorRadius]).^2);
if ~isSphere && opt.enhancingSourceShape == "spheroid" && nearestProxyLevel >= 1
    error('lap3d_build_two_body_grids:ProxySegmentNotEnclosed', ...
        'The stretched spheroid does not enclose the nearest proxy point.');
end

clearance = ellipsoid_surface_clearance(E,parentPoints);
chopDepth = opt.sourceOffset+opt.sourceChopOffset;
keep = clearance < chopDepth;
cap = parentPoints(keep,:);
if isempty(cap)
    error('lap3d_build_two_body_grids:EmptySourceCap', ...
        ['Chopping removed every enhancing source node. ', ...
         'Increase sourceChopOffset, sourceBeta, or sourceSphereRadius.']);
end

metadata.shape = opt.enhancingSourceShape;
metadata.distribution = opt.enhancingSourceDistribution;
metadata.requestedResolution = opt.fineSourceResolution;
metadata.requestedDesignCount = opt.enhancingSourceCount;
metadata.parentNodeCount = size(parentPoints,1);
metadata.chopOffset = opt.sourceChopOffset;
metadata.chopDepth = chopDepth;
metadata.focalHalfLength = focalHalfLength;
metadata.widthScale = opt.sourceWidthScale;
metadata.parentFocalHalfLength = sqrt(max(0,majorRadius^2-minorRadius^2));
metadata.parentPoints = parentPoints;
metadata.parentAxes = [minorRadius;minorRadius;majorRadius];
metadata.parentCentre = proxyCentre;
metadata.parentRotation = proxyRotation;
metadata.nearFocus = nearFocus;
metadata.innerFocus = innerFocus;
metadata.nearAccumulationPoint = nearAccumulationPoint;
metadata.radialSurfacePoint = radialSurfacePoint;
metadata.distanceToSurface = distanceToSurface;
metadata.outerTip = proxyCentre+proxyRotation(:,3)*majorRadius;
metadata.nearestBasicProxyPoint = nearestBasicProxyPoint;
metadata.proxyReferencePoint = proxyReferencePoint;
metadata.proxyReferenceIsSampled = proxyReferenceIsSampled;
metadata.proxyRadius = proxyRadius;
metadata.isMirrored = false;
metadata.requestedTipExtension = requestedTipExtension;
metadata.requiredMajorRadius = requiredMajorRadius;
metadata.effectiveBeta = tipExtension/distanceToSurface;
metadata.nearestProxyLevel = nearestProxyLevel;
metadata.maxAccumulationLevel = max(accumulationLevel);
metadata.maxPhysicalLevel = max(physicalLevel);
metadata.retainedCountBeforeDeduplication = nnz(keep);
metadata.tipExtension = tipExtension;
metadata.minimumRetainedClearance = min(clearance(keep));
metadata.maximumRetainedClearance = max(clearance(keep));
end


function Q = frame_with_z_axis(zAxis,trialDirection)
%FRAME_WITH_Z_AXIS constructs a proper orthogonal frame for a prolate proxy.

zAxis = zAxis/norm(zAxis);
xAxis = trialDirection-zAxis*(zAxis.'*trialDirection);
if norm(xAxis) < 1e-12
    trialDirection = [0;0;1];
    if abs(trialDirection.'*zAxis) > 0.9
        trialDirection = [0;1;0];
    end
    xAxis = trialDirection-zAxis*(zAxis.'*trialDirection);
end
xAxis = xAxis/norm(xAxis);
yAxis = cross(zAxis,xAxis);
yAxis = yAxis/norm(yAxis);
Q = [xAxis,yAxis,zAxis];
end


function clearance = ellipsoid_surface_clearance(E,points)
%ELLIPSOID_SURFACE_CLEARANCE measures interior points to the physical surface.

localPoints = E.Q.'*(points.'-E.c);
if max(E.axes)-min(E.axes) < 1e-12*max(E.axes)
    clearance = E.axes(1)-vecnorm(localPoints,2,1).';
    return
end
clearance = zeros(size(points,1),1);
for k = 1:size(points,1)
    [~,nearestPoint] = dist_ellipsoid(localPoints(:,k),E.axes);
    clearance(k) = norm(localPoints(:,k)-nearestPoint);
end
end


function points = remove_near_duplicates(candidates,existing,tolerance)
%REMOVE_NEAR_DUPLICATES removes repeated candidates without extra toolboxes.

points = zeros(0,3);
accepted = existing;
for k = 1:size(candidates,1)
    if isempty(accepted)
        separated = true;
    else
        difference = accepted-candidates(k,:);
        separated = min(sum(difference.^2,2)) > tolerance^2;
    end
    if separated
        points(end+1,:) = candidates(k,:); %#ok<AGROW>
        accepted(end+1,:) = candidates(k,:); %#ok<AGROW>
    end
end
end


function level = ellipsoid_level(E,points)
%ELLIPSOID_LEVEL evaluates the physical ellipsoid's implicit level function.

localPoints = E.Q.'*(points.'-E.c);
scaledPoints = localPoints./E.axes;
level = sum(scaledPoints.^2,1).';
end


function validate_grid_points(grid,bodies)
%VALIDATE_GRID_POINTS checks point ownership and source/target placement.

if any(~isfinite(grid.sources),'all') || any(~isfinite(grid.collocation),'all')
    error('lap3d_build_two_body_grids:NonfiniteGrid', ...
        'All source and collocation coordinates must be finite.');
end

for body = 1:2
    sourceLevel = ellipsoid_level(bodies{body}, ...
        grid.sources(grid.sourceBody == body,:));
    collocationLevel = ellipsoid_level(bodies{body}, ...
        grid.collocation(grid.collocationBody == body,:));
    if max(sourceLevel) >= 1
        error('lap3d_build_two_body_grids:SourceOutsideBody', ...
            'Body %d has a source outside or on the physical surface.',body);
    end
    if max(abs(collocationLevel-1)) > 5e-9
        error('lap3d_build_two_body_grids:TargetOffSurface', ...
            'Body %d has a collocation point off the physical surface.',body);
    end
end
end
