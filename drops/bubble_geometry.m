function geometry = bubble_geometry(opt)
%BUBBLE_GEOMETRY Interior Stokeslets and off-edge mixed-boundary nodes.
if nargin == 0
    bubble_self_test();
    geometry = [];
    return
end
assert(isfield(opt,'geometryOptions'),'Complete options with bubble_options.');
if opt.n_level > 0
    assert(opt.phi>0 && opt.phi<pi, ...
        'Cone rings require a nondegenerate cap edge.');
    % Hold source sites and ring offsets fixed across basis choices. Add
    % collocation azimuths for each extra cone block before building rings.
    baseOptions = opt.geometryOptions;
    baseOptions.useExtraCollocation = false;
    base = stagnant_cap_geometry(baseOptions);
    coneCount = opt.n_level*opt.n_az_source;
    familyCount = numel(opt.refinementFamilies)+ ...
        double(opt.derivative_direction~="none");
    ringCount = opt.ringLevels;
    baselinePoints = ceil(opt.collocationRatio* ...
        (base.N+coneCount));
    baselineAzimuths = max([opt.n_az_colloc,4, ...
        ceil((baselinePoints-base.collocationCount)/(2*ringCount))]);
    addedAzimuths = ceil(opt.collocationRatio* ...
        (familyCount-1)*coneCount/(2*ringCount));
    requiredPoints = ceil(opt.collocationRatio* ...
        (base.N+familyCount*coneCount));
    neededAzimuths = ceil((requiredPoints-base.collocationCount)/ ...
        (2*ringCount));
    geometryOptions = opt.geometryOptions;
    geometryOptions.n_az_colloc = max( ...
        baselineAzimuths+addedAzimuths,neededAzimuths);
    raw = stagnant_cap_refined_geometry(geometryOptions);
    sourceRingIndex = raw.innerSourceRingIndex;
    collocationRingIndex = raw.collocationRingIndex;
    coneCount = size(raw.innerAddedSources,1);
    collocationRingCount = raw.collocationRingCount;
    assert(isequal(collocationRingCount,[opt.ringLevels opt.ringLevels]));
else
    raw = stagnant_cap_geometry(opt.geometryOptions);
    sourceRingIndex = zeros(raw.N,1);
    coneCount = 0;
    collocationRingIndex = zeros(raw.collocationCount,1);
    % Reposition only the added rings; preserve the uniform base design.
    nRingLevel = numel(raw.ringOffsets);
    if nRingLevel>0
        offsets = opt.edgeBandAngle* ...
            opt.collocClusterRatio.^(0:nRingLevel-1);
        nAz = raw.ringAzimuthCount;
        for level = 1:nRingLevel
            for side = [-1 1]
                first = raw.baseCount+(level-1)*2*nAz+ ...
                    (side+1)*nAz/2+1;
                rows = first:first+nAz-1;
                theta = opt.phi+side*offsets(level);
                phase = mod(level/2+(side+1)/4,1);
                psi = 2*pi*((0:nAz-1)+phase)/nAz;
                unit = [sin(theta)*cos(psi(:)), ...
                    sin(theta)*sin(psi(:)),cos(theta)*ones(nAz,1)];
                raw.points(rows,:) = opt.centre+opt.radius*unit;
                raw.normals(rows,:) = unit;
                raw.theta(rows) = theta;
                raw.cap(rows) = side<0;
            end
        end
        raw.mobile = ~raw.cap;
        raw.edgeDistance = opt.radius*abs(raw.theta-opt.phi);
        raw.ringOffsets = offsets(:);
    end
    for level = 1:nRingLevel
        first = raw.baseCount+(level-1)*2*raw.ringAzimuthCount+1;
        last = first+2*raw.ringAzimuthCount-1;
        collocationRingIndex(first:last) = level;
    end
    collocationRingCount = [nRingLevel,nRingLevel];
end
geometry.sourcePoints = raw.innerSources;
geometry.sourceCount = size(geometry.sourcePoints,1);
geometry.baseSourceCount = raw.N;
geometry.coneSourceCount = coneCount;
if opt.n_level>0
    geometry.coneFamilyCount = familyCount;
else
    geometry.coneFamilyCount = 0;
end
geometry.basisPointCount = geometry.baseSourceCount+ ...
    geometry.coneFamilyCount*geometry.coneSourceCount;
geometry.coefficientCount = 3*geometry.basisPointCount;
geometry.sourceRingCount = opt.n_level;
geometry.sourceRingIndex = sourceRingIndex;
geometry.sourceRegion = repmat("base",geometry.sourceCount,1);
geometry.sourceRegion(sourceRingIndex>0) = "cone";
geometry.points = raw.points;
geometry.normals = raw.normals;
geometry.theta = raw.theta;
geometry.cap = raw.cap;
geometry.mobile = raw.mobile;
geometry.collocationCount = raw.collocationCount;
geometry.capCount = nnz(raw.cap);
geometry.mobileCount = nnz(raw.mobile);
geometry.baseCollocationCount = raw.baseCount;
geometry.collocationRingIndex = collocationRingIndex;
geometry.collocationRingCount = collocationRingCount;
geometry.collocationRegion = repmat("mobile",raw.collocationCount,1);
geometry.collocationRegion(raw.cap) = "cap";
geometry.ringOffsets = raw.ringOffsets;
geometry.ringAzimuthCount = raw.ringAzimuthCount;
geometry.edgeDistance = raw.edgeDistance;
geometry.minimumSourceDistance = min(opt.radius- ...
    vecnorm(geometry.sourcePoints-opt.centre,2,2));
geometry.minimumEdgeDistance = min(raw.edgeDistance);
geometry.minimumAddedEdgeDistance = NaN;
if ~isempty(raw.ringOffsets)
    geometry.minimumAddedEdgeDistance = opt.radius*min(raw.ringOffsets);
end
assert(geometry.sourceCount==geometry.baseSourceCount+ ...
    opt.n_level*opt.n_az_source);
assert(coneCount==opt.n_level*opt.n_az_source);
assert(geometry.capCount+geometry.mobileCount==geometry.collocationCount);
assert(geometry.minimumSourceDistance>0);
assert(geometry.collocationCount>geometry.sourceCount, ...
    'Need more boundary points than Stokeslet points.');
if opt.n_level>0
    assert(3*geometry.collocationCount>= ...
        opt.collocationRatio*geometry.coefficientCount);
end
assert(all(geometry.edgeDistance>1e-12*opt.radius) || ...
    opt.phi==0 || opt.phi==pi);
end

function bubble_self_test()
% Verify five interior cone rings, counts, uniqueness, and both geometry views.
o = bubble_options(struct('sourceRequest',50,'phi',pi/4, ...
    'n_level',5,'n_az_source',8,'n_az_colloc',12));
g = bubble_geometry(o);
assert(g.coneSourceCount==40 && g.sourceRingCount==5);
assert(isequal(g.collocationRingCount,[10 10]));
assert(nnz(g.sourceRingIndex==5)==8);
assert(all(vecnorm(g.sourcePoints-o.centre,2,2)<o.radius));
mixedOptions = o;
mixedOptions.refinementFamilies = ["S","T_t","T_n","T_n2"];
mixedOptions = bubble_options(mixedOptions);
mixed = bubble_geometry(mixedOptions);
assert(isequal(g.sourcePoints,mixed.sourcePoints) && ...
    isequal(g.ringOffsets,mixed.ringOffsets) && ...
    isequal(g.collocationRingCount,mixed.collocationRingCount));
assert(mixed.ringAzimuthCount>g.ringAzimuthCount && ...
    mixed.coefficientCount==3*(g.baseSourceCount+4*g.coneSourceCount));
assert(mixed.collocationCount>=g.collocationCount+ ...
    ceil(o.collocationRatio*3*g.coneSourceCount));
for j = 1:g.sourceCount-1
    assert(all(vecnorm(g.sourcePoints(j+1:end,:)- ...
        g.sourcePoints(j,:),2,2)>1e-11*o.radius));
end
fprintf('Bubble geometry: %d base + %d cone sources; %d cap + %d mobile nodes.\n', ...
    g.baseSourceCount,g.coneSourceCount,g.capCount,g.mobileCount);
off = o;
off.n_level = 0;
off = bubble_options(off);
uniform = bubble_geometry(off);
assert(uniform.coneSourceCount==0 && ...
    isequal(uniform.sourcePoints,g.sourcePoints(1:g.baseSourceCount,:)));
uniformDefault = bubble_options(struct('sourceRequest',50,'phi',pi/4, ...
    'n_level',0,'ringLevels',4));
uniformBaseline = bubble_geometry(uniformDefault);
assert(abs(uniformBaseline.minimumAddedEdgeDistance/ ...
    uniformDefault.radius-uniformDefault.edgeBandAngle*0.5^3)<1e-13);
uniformControl = bubble_options(struct('sourceRequest',50,'phi',pi/4, ...
    'n_level',0,'ringLevels',4,'collocClusterRatio',0.3));
uniformClustered = bubble_geometry(uniformControl);
assert(abs(uniformClustered.minimumAddedEdgeDistance/ ...
    uniformControl.radius-uniformControl.edgeBandAngle*0.3^3)<1e-13);
assert(isequal(uniformBaseline.sourcePoints,uniformClustered.sourcePoints));
assert(isequal(uniformBaseline.points(1:uniformBaseline.baseCollocationCount,:), ...
    uniformClustered.points(1:uniformClustered.baseCollocationCount,:)));
[f3,f2] = bubble_plot_geometry(g,o,'on');
assert(isgraphics(f3) && isgraphics(f2));
fprintf('bubble_geometry self-test passed.\n');
end
