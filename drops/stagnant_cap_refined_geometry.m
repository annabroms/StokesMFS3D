function geometry = stagnant_cap_refined_geometry(opt)
%STAGNANT_CAP_REFINED_GEOMETRY Add conical sources and graded surface rings.
% Base nodes retain their original order. Ring indices are zero for base
% nodes and 1:n_level (sources) or 1:n_colloc_level (surface nodes).
if nargin==0
    self_test();
    return
end
opt = stagnant_cap_options(opt);
assert(opt.phi>0 && opt.phi<pi, ...
    'Source refinement requires a nondegenerate cap edge.');
assert(opt.collocAngularExtent>0 && ...
    opt.collocAngularExtent<min(opt.phi,pi-opt.phi), ...
    'Collocation refinement must stay within both surface regions.');
baseOpt = opt;
baseOpt.useExtraCollocation = false;
geometry = stagnant_cap_geometry(baseOpt);
a = opt.radius;
c = opt.centre;
nLevel = opt.n_level;
nAzSource = opt.n_az_source;
nFamily = numel(opt.refinementFamilies)+double(opt.derivative_direction~="none");
nCollocLevel = ceil(opt.colloc_factor*nLevel*nFamily);
innerAdded = zeros(nLevel*nAzSource,3);
outerAdded = zeros(nLevel*nAzSource,3);
innerDistance = opt.innerSourceStartDistance* ...
    opt.sourceClusterRatio.^(0:nLevel-1);
outerDistance = opt.outerSourceStartDistance* ...
    opt.sourceClusterRatio.^(0:nLevel-1);
tol = 1e-11*a;
for k = 1:nLevel
    rows = (k-1)*nAzSource+(1:nAzSource);
    innerAdded(rows,:) = place_ring(c,a-innerDistance(k), ...
        opt.phi,nAzSource,[geometry.innerSources;innerAdded(1:rows(1)-1,:)],tol);
    outerAdded(rows,:) = place_ring(c,a+outerDistance(k), ...
        opt.phi,nAzSource,[geometry.outerSources;outerAdded(1:rows(1)-1,:)],tol);
end
baseN = geometry.N;
geometry.innerSources = [geometry.innerSources;innerAdded];
geometry.outerSources = [geometry.outerSources;outerAdded];
geometry.Nin = size(geometry.innerSources,1);
geometry.Nex = size(geometry.outerSources,1);
geometry.innerSourceRingIndex = [zeros(baseN,1); ...
    reshape(repelem((1:nLevel).',nAzSource),[],1)];
geometry.outerSourceRingIndex = geometry.innerSourceRingIndex;
geometry.sourceRegion = [repmat("inside",geometry.Nin,1); ...
    repmat("outside",geometry.Nex,1)];
geometry.sourceRingIndex = [geometry.innerSourceRingIndex; ...
    geometry.outerSourceRingIndex];
geometry.innerBaseSources = geometry.innerSources(1:baseN,:);
geometry.outerBaseSources = geometry.outerSources(1:baseN,:);
geometry.innerAddedSources = innerAdded;
geometry.outerAddedSources = outerAdded;
geometry.sourceDistances.inner = innerDistance(:);
geometry.sourceDistances.outer = outerDistance(:);
geometry.innerDerivativeDistances = reshape( ...
    repelem(innerDistance(:),nAzSource),[],1);
geometry.outerDerivativeDistances = reshape( ...
    repelem(outerDistance(:),nAzSource),[],1);
geometry.innerDerivativeDirections = derivative_directions( ...
    innerAdded,c,opt.phi,opt.derivative_direction);
geometry.outerDerivativeDirections = derivative_directions( ...
    outerAdded,c,opt.phi,opt.derivative_direction);
geometry.enhancingFamilyCount = nFamily;
geometry.sourceRingCount = [nLevel,nLevel];
geometry.addedSourceCount = [size(innerAdded,1),size(outerAdded,1)];

baseBasis = struct('family',"S",'points',zeros(0,3), ...
    'normalDirections',zeros(0,3),'scale',1, ...
    'derivativeDistances',zeros(0,1));
geometry.exteriorBasis = baseBasis;
geometry.exteriorBasis.points = geometry.innerBaseSources;
geometry.interiorBasis = baseBasis;
geometry.interiorBasis.points = geometry.outerBaseSources;
for family = opt.refinementFamilies
    block = baseBasis;
    block.family = family;
    block.points = innerAdded;
    block.normalDirections = source_directions( ...
        innerAdded,c,family,opt.use_orth_dir_images);
    if family=="Q_s", block.scale = a^2; end
    if ~isempty(innerAdded)
        geometry.exteriorBasis(end+1,1) = block;
    end
    block.points = outerAdded;
    block.normalDirections = source_directions( ...
        outerAdded,c,family,opt.use_orth_dir_images);
    if ~isempty(outerAdded)
        geometry.interiorBasis(end+1,1) = block;
    end
end
if opt.derivative_direction~="none" && ~isempty(innerAdded)
    block = baseBasis;
    block.family = "S_y";
    block.points = innerAdded;
    block.normalDirections = geometry.innerDerivativeDirections;
    block.derivativeDistances = geometry.innerDerivativeDistances;
    geometry.exteriorBasis(end+1,1) = block;
    block.points = outerAdded;
    block.normalDirections = geometry.outerDerivativeDirections;
    block.derivativeDistances = geometry.outerDerivativeDistances;
    geometry.interiorBasis(end+1,1) = block;
end
geometry.exteriorCoefficientCount = 3*sum(arrayfun( ...
    @(b) size(b.points,1),geometry.exteriorBasis));
geometry.interiorCoefficientCount = 3*sum(arrayfun( ...
    @(b) size(b.points,1),geometry.interiorBasis));

neededM = ceil(opt.collocationRatio* ...
    (geometry.exteriorCoefficientCount+geometry.interiorCoefficientCount)/6);
if nCollocLevel>0
    nAzColloc = max(opt.n_az_colloc, ...
        ceil((neededM-geometry.baseCount)/(2*nCollocLevel)));
    nAzColloc = max(nAzColloc,4);
else
    nAzColloc = 0;
end

offsets = opt.collocAngularExtent* ...
    opt.collocClusterRatio.^(0:nCollocLevel-1);
ringIndex = zeros(geometry.baseCount+2*nCollocLevel*nAzColloc,1);
unit = (geometry.points-c)/a;
for k = 1:nCollocLevel
    for side = [-1,1]
        if side<0
            theta = opt.phi-offsets(k);
        else
            theta = opt.phi+offsets(k);
        end
        x = place_ring(c,a,theta,nAzColloc,c+a*unit,tol);
        unit = [unit;(x-c)/a]; %#ok<AGROW>
        rows = size(unit,1)-nAzColloc+(1:nAzColloc);
        ringIndex(rows) = k;
    end
end
geometry.points = c+a*unit;
geometry.normals = unit;
geometry.theta = acos(max(-1,min(1,unit(:,3))));
geometry.cap = geometry.theta<opt.phi;
geometry.mobile = ~geometry.cap;
geometry.collocationRegion = repmat("mobile",size(unit,1),1);
geometry.collocationRegion(geometry.cap) = "cap";
geometry.edgeDistance = a*abs(geometry.theta-opt.phi);
geometry.collocationRingIndex = ringIndex;
geometry.baseCollocation = ringIndex==0;
geometry.ringCollocation = ringIndex>0;
geometry.collocationCount = size(unit,1);
geometry.baseCollocationPoints = geometry.points(geometry.baseCollocation,:);
geometry.addedCollocationPoints = geometry.points(geometry.ringCollocation,:);
geometry.ringOffsets = offsets(:);
geometry.ringAzimuthCount = nAzColloc;
geometry.collocationRingCount = [nCollocLevel,nCollocLevel];
geometry.addedCollocationCount = [nnz(geometry.ringCollocation & geometry.cap), ...
    nnz(geometry.ringCollocation & geometry.mobile)];
assert(all(vecnorm(geometry.innerSources-c,2,2)<a) && ...
    all(vecnorm(geometry.outerSources-c,2,2)>a));
assert(all(geometry.edgeDistance>tol));
assert(6*geometry.collocationCount>=opt.collocationRatio* ...
    (geometry.exteriorCoefficientCount+geometry.interiorCoefficientCount));
assert_distinct(geometry.innerSources,tol);
assert_distinct(geometry.outerSources,tol);
assert_distinct(geometry.points,tol);
end

function directions = derivative_directions(points,c,phi,choice)
%DERIVATIVE_DIRECTIONS Use the edge frame with one sign on both source sides.
directions = zeros(size(points));
if isempty(points) || choice=="none", return; end
psi = atan2(points(:,2)-c(2),points(:,1)-c(1));
if choice=="radial"
    directions = [sin(phi)*cos(psi),sin(phi)*sin(psi), ...
        cos(phi)*ones(size(psi))];
else
    directions = [cos(phi)*cos(psi),cos(phi)*sin(psi), ...
        -sin(phi)*ones(size(psi))];
end
end

function ring = place_ring(c,r,theta,nAz,existing,tol)
%PLACE_RING Rotate a full latitude ring until it avoids existing nodes.
for attempt = 0:16
    phase = mod(0.173+attempt*sqrt(2)/17,1);
    psi = 2*pi*((0:nAz-1)+phase)/nAz;
    ring = c+r*[sin(theta)*cos(psi(:)), ...
        sin(theta)*sin(psi(:)),cos(theta)*ones(nAz,1)];
    if isempty(existing), return; end
    distinct = true;
    for j = 1:nAz
        if any(sum((existing-ring(j,:)).^2,2)<=tol^2)
            distinct = false;
            break
        end
    end
    if distinct, return; end
end
error('stagnant_cap_refined_geometry:CoincidentNodes', ...
    'A source or collocation ring could not avoid an existing node.');
end

function directions = source_directions(points,c,family,useOrthogonal)
%SOURCE_DIRECTIONS Reuse image-source orientation rules on radial rays.
directions = zeros(0,3);
if isempty(points), return; end
ray = points-c;
ray = ray./vecnorm(ray,2,2);
if ismember(family,["T","T_n","T_n2","T_t"])
    directions = image_stresslet_directions( ...
        ray,useOrthogonal,family);
elseif family=="Q_s"
    directions = ray;
end
end

function assert_distinct(points,tol)
%ASSERT_DISTINCT Reject coincident three-dimensional nodes in one geometry.
for j = 1:size(points,1)-1
    assert(~any(sum((points(j+1:end,:)-points(j,:)).^2,2)<=tol^2), ...
        'Merged geometry contains duplicate points.');
end
end

function self_test()
%SELF_TEST Verify five conical levels, counts, and both requested views.
opt = stagnant_cap_options(struct('sourceRequest',50,'phi',pi/4, ...
    'n_level',5,'n_az_source',64,'n_az_colloc',64));
g = stagnant_cap_refined_geometry(opt);
assert(isequal(g.sourceRingCount,[5 5]) && ...
    isequal(g.collocationRingCount,[10 10]));
assert(isequal(g.addedSourceCount,[320 320]) && ...
    isequal(g.addedCollocationCount,[640 640]));
assert(abs(g.sourceDistances.inner(end)-0.108045*opt.radius)<1e-12);
assert(all(g.innerSourceRingIndex(1:g.N)==0));
assert(all(g.outerSourceRingIndex(1:g.N)==0));
assert(all(g.collocationRingIndex(1:g.baseCount)==0));
assert(max(abs(vecnorm(g.innerAddedSources-opt.centre,2,2)- ...
    repelem((opt.radius-g.sourceDistances.inner),opt.n_az_source)))<1e-12);
firstRing = place_ring(opt.centre,opt.radius,opt.phi,8,[],1e-11);
rotated = place_ring(opt.centre,opt.radius,opt.phi,8, ...
    firstRing(1,:),1e-11);
assert(min(vecnorm(rotated-firstRing(1,:),2,2))>1e-11);
off = opt;
off.n_level = 0;
noRings = stagnant_cap_refined_geometry(off);
assert(isempty(noRings.innerAddedSources) && ...
    isempty(noRings.addedCollocationPoints) && ...
    noRings.collocationCount==noRings.baseCount);
multi = opt;
multi.n_level = 1;
multi.n_az_source = 6;
multi.n_az_colloc = 4;
multi.refinementFamilies = "SDRT_nT_tT_n2Q_s";
mixed = stagnant_cap_refined_geometry(multi);
assert(mixed.ringAzimuthCount>=multi.n_az_colloc && ...
    isequal(mixed.collocationRingCount,[14 14]) && ...
    numel(mixed.exteriorBasis)==1+7);
directions = {mixed.exteriorBasis(5).normalDirections, ...
    mixed.exteriorBasis(6).normalDirections, ...
    mixed.exteriorBasis(7).normalDirections};
for j = 1:3
    for k = j+1:3
        assert(max(abs(sum(directions{j}.*directions{k},2)))<1e-12);
    end
end
assert(max(abs(vecnorm(mixed.exteriorBasis(end).normalDirections, ...
    2,2)-1))<1e-12);
radialOpt = opt;
radialOpt.derivative_direction = "radial";
radial = stagnant_cap_refined_geometry(radialOpt);
coneOpt = opt;
coneOpt.derivative_direction = "cone_normal";
cone = stagnant_cap_refined_geometry(coneOpt);
assert(isequal(radial.innerAddedSources,cone.innerAddedSources) && ...
    isequal(radial.outerAddedSources,cone.outerAddedSources) && ...
    isequal(radial.points,cone.points) && ...
    isequal(radial.collocationRingCount,[20 20]));
for side = ["inner","outer"]
    qr = radial.(side+"DerivativeDirections");
    qt = cone.(side+"DerivativeDirections");
    assert(max(abs(vecnorm(qr,2,2)-1))<1e-12 && ...
        max(abs(vecnorm(qt,2,2)-1))<1e-12 && ...
        max(abs(sum(qr.*qt,2)))<1e-12);
end
psi = atan2(radial.innerAddedSources(:,2)-opt.centre(2), ...
    radial.innerAddedSources(:,1)-opt.centre(1));
expectedRadial = [sin(opt.phi)*cos(psi),sin(opt.phi)*sin(psi), ...
    cos(opt.phi)*ones(size(psi))];
expectedCone = [cos(opt.phi)*cos(psi),cos(opt.phi)*sin(psi), ...
    -sin(opt.phi)*ones(size(psi))];
assert(max(abs(radial.innerDerivativeDirections-expectedRadial),[],'all') ...
    <1e-12);
assert(max(abs(cone.innerDerivativeDirections-expectedCone),[],'all') ...
    <1e-12);
assert(isequal(radial.innerDerivativeDistances,reshape(repelem( ...
    radial.sourceDistances.inner,opt.n_az_source),[],1)));
assert(max(abs(radial.innerDerivativeDirections- ...
    radial.outerDerivativeDirections),[],'all')<1e-12);
assert(max(abs(cone.innerDerivativeDirections- ...
    cone.outerDerivativeDirections),[],'all')<1e-12);

f3 = figure('Name','Refined stagnant-cap geometry: 3D','Color','w');
ax = axes(f3); hold(ax,'on');
[sx,sy,sz] = sphere(48);
surf(ax,opt.centre(1)+opt.radius*sx,opt.centre(2)+opt.radius*sy, ...
    opt.centre(3)+opt.radius*sz,'FaceAlpha',0.12,'EdgeColor','none');
psi = linspace(0,2*pi,300);
edge = opt.centre+opt.radius*[sin(opt.phi)*cos(psi(:)), ...
    sin(opt.phi)*sin(psi(:)),cos(opt.phi)*ones(numel(psi),1)];
plot3(ax,edge(:,1),edge(:,2),edge(:,3),'k-','LineWidth',2);
sets = {radial.innerBaseSources,radial.outerBaseSources, ...
    radial.innerAddedSources,radial.outerAddedSources, ...
    radial.points(radial.ringCollocation & radial.cap,:), ...
    radial.points(radial.ringCollocation & radial.mobile,:)};
colors = {[0.1 0.3 0.8],[0.8 0.35 0.1], ...
    [0 0.6 0.9],[0.9 0.15 0.1], ...
    [0.8 0.1 0.6],[0.1 0.7 0.3]};
for j = 1:numel(sets)
    p = sets{j};
    scatter3(ax,p(:,1),p(:,2),p(:,3),9,colors{j},'filled');
end
arrowIndices = 1:8:opt.n_az_source;
arrowScale = 0.22*opt.radius;
for side = 1:2
    if side==1, points = radial.innerAddedSources; ...
    else, points = radial.outerAddedSources; end
    p = points(arrowIndices,:);
    qr = radial.innerDerivativeDirections(arrowIndices,:);
    qt = cone.innerDerivativeDirections(arrowIndices,:);
    quiver3(ax,p(:,1),p(:,2),p(:,3),arrowScale*qr(:,1), ...
        arrowScale*qr(:,2),arrowScale*qr(:,3),0,'b','LineWidth',1.2);
    quiver3(ax,p(:,1),p(:,2),p(:,3),arrowScale*qt(:,1), ...
        arrowScale*qt(:,2),arrowScale*qt(:,3),0,'m','LineWidth',1.2);
end
axis(ax,'equal'); view(ax,35,25); grid(ax,'on');
legend(ax,{'interface','cap edge','base inside','base outside', ...
    'added inside','added outside','extra cap','extra mobile', ...
    'radial derivative','cone-normal derivative'});

fm = figure('Name','Refined stagnant-cap geometry: meridian','Color','w');
ax = axes(fm); hold(ax,'on');
theta = linspace(0,pi,300);
plot(ax,opt.radius*sin(theta),opt.radius*cos(theta),'k-');
plot(ax,opt.radius*sin(opt.phi),opt.radius*cos(opt.phi), ...
    'kp','MarkerFaceColor','k');
for j = 1:numel(sets)
    p = sets{j}-opt.centre;
    scatter(ax,hypot(p(:,1),p(:,2)),p(:,3),12,colors{j},'filled');
end
for side = 1:2
    if side==1, points = radial.innerAddedSources; ...
    else, points = radial.outerAddedSources; end
    p = points(arrowIndices,:)-opt.centre;
    psi = atan2(p(:,2),p(:,1));
    qr = radial.innerDerivativeDirections(arrowIndices,:);
    qt = cone.innerDerivativeDirections(arrowIndices,:);
    qrr = qr(:,1).*cos(psi)+qr(:,2).*sin(psi);
    qtr = qt(:,1).*cos(psi)+qt(:,2).*sin(psi);
    quiver(ax,hypot(p(:,1),p(:,2)),p(:,3),arrowScale*qrr, ...
        arrowScale*qr(:,3),0,'b','LineWidth',1.2);
    quiver(ax,hypot(p(:,1),p(:,2)),p(:,3),arrowScale*qtr, ...
        arrowScale*qt(:,3),0,'m','LineWidth',1.2);
end
axis(ax,'equal'); grid(ax,'on');
xlabel(ax,'distance from axis'); ylabel(ax,'z minus centre');
legend(ax,{'interface','cap edge','base inside','base outside', ...
    'added inside','added outside','extra cap','extra mobile', ...
    'radial derivative','cone-normal derivative'});
fprintf(['Refined geometry: source rings %d inside + %d outside, ', ...
    'points %d + %d = %d; collocation rings %d cap + %d mobile, ', ...
    'points %d + %d = %d.\n'],g.sourceRingCount, ...
    g.addedSourceCount,sum(g.addedSourceCount), ...
    g.collocationRingCount,g.addedCollocationCount, ...
    sum(g.addedCollocationCount));
fprintf('Derivative geometry: %d collocation rings per region.\n', ...
    radial.collocationRingCount(1));
end
