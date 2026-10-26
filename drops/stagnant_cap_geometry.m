function geometry = stagnant_cap_geometry(opt)
%STAGNANT_CAP_GEOMETRY Uniform proxy designs and graded surface collocation.
if nargin == 0
    self_test();
    return
end
opt = stagnant_cap_options(opt);
[sourceUnit,~] = stagnant_cap_design(opt.sourceRequest);
N = size(sourceUnit,1);
geometry.innerSources = opt.centre + opt.radius*opt.innerProxyRatio*sourceUnit;
geometry.outerSources = opt.centre + opt.radius*opt.outerProxyRatio*sourceUnit;
geometry.N = N; % Number of Stokeslet sources on each proxy sphere.
[base,~] = stagnant_cap_design(ceil(opt.collocationRatio*N));
baseTheta = acos(max(-1,min(1,base(:,3))));
hasPhysicalEdge = opt.phi>0 && opt.phi<pi;
if hasPhysicalEdge
    % Neither mode may collocate exactly where the boundary law changes.
    keepOffEdge = abs(baseTheta-opt.phi)>1e-12;
    base = base(keepOffEdge,:);
end
hasEdge = hasPhysicalEdge && opt.useExtraCollocation && ...
    opt.ringLevels>0 && opt.edgeBandAngle>0;
if hasEdge
    % Preserve every off-edge base node; rings only add edge resolution.
    offsets = opt.edgeBandAngle*2.^(-(0:opt.ringLevels-1));
    if isempty(opt.ringAzimuthCount)
        spacing = sqrt(4*pi/ceil(opt.collocationRatio*N));
        nAz = max(12,ceil(2*pi*sin(opt.phi)/spacing));
    else
        nAz = opt.ringAzimuthCount;
    end
    rings = zeros(2*opt.ringLevels*nAz,3);
    row = 0;
    for level = 1:opt.ringLevels
        for side = [-1 1]
            theta = opt.phi+side*offsets(level);
            phase = mod(level/2+(side+1)/4,1);
            psi = 2*pi*((0:nAz-1)+phase)/nAz;
            selected = row+(1:nAz);
            rings(selected,:) = [sin(theta)*cos(psi(:)), ...
                sin(theta)*sin(psi(:)),cos(theta)*ones(nAz,1)];
            row = row+nAz;
        end
    end
    unit = [base;rings];
else
    offsets = [];
    nAz = 0;
    unit = base;
end
geometry.points = opt.centre+opt.radius*unit;
geometry.normals = unit;
geometry.theta = acos(max(-1,min(1,unit(:,3))));
geometry.cap = geometry.theta<opt.phi;
geometry.mobile = ~geometry.cap;
geometry.edgeDistance = opt.radius*abs(geometry.theta-opt.phi);
geometry.bandAngle = opt.edgeBandAngle;
geometry.ringOffsets = offsets(:);
geometry.ringAzimuthCount = nAz;
geometry.baseCount = size(base,1);
geometry.collocationCount = size(unit,1);
geometry.baseCollocation = false(geometry.collocationCount,1);
geometry.baseCollocation(1:geometry.baseCount) = true;
geometry.ringCollocation = ~geometry.baseCollocation;
assert(geometry.collocationCount>N);
assert(all(geometry.edgeDistance>1e-12*opt.radius) || ~hasPhysicalEdge);
assert(all(vecnorm(geometry.innerSources-opt.centre,2,2)<opt.radius));
assert(all(vecnorm(geometry.outerSources-opt.centre,2,2)>opt.radius));
end

function self_test()
%SELF_TEST Verify uniform proxies, both edge sides, and excluded edge nodes.
o = stagnant_cap_options(struct('sourceRequest',50,'ringLevels',3));
g = stagnant_cap_geometry(o);
assert(g.N>=50 && g.collocationCount>g.N);
assert(nnz(g.cap)>0 && nnz(g.mobile)>0);
assert(numel(g.ringOffsets)==3);
assert(g.collocationCount==g.baseCount+2*o.ringLevels*g.ringAzimuthCount);
assert(nnz(g.ringCollocation)==2*o.ringLevels*g.ringAzimuthCount);
ringThetas = [o.phi-g.ringOffsets;o.phi+g.ringOffsets];
assert(max(abs(sort(g.theta(g.ringCollocation))- ...
    sort(repelem(ringThetas(:),g.ringAzimuthCount))))<1e-12);
assert(max(abs(vecnorm(g.normals,2,2)-1))<1e-12);
assert(max(abs(vecnorm(g.innerSources-o.centre,2,2)- ...
    o.radius*o.innerProxyRatio))<1e-12);
o.useExtraCollocation = false;
off = stagnant_cap_geometry(o);
assert(isempty(off.ringOffsets) && off.ringAzimuthCount==0);
assert(off.collocationCount==off.baseCount);
assert(off.collocationCount>off.N && all(off.edgeDistance>1e-12*o.radius));
assert(isequal(g.points(1:g.baseCount,:),off.points));
assert(isequal(g.cap(1:g.baseCount),off.cap));
assert(~any(off.ringCollocation) && all(off.baseCollocation));
assert(isequal(off.innerSources,g.innerSources) && ...
    isequal(off.outerSources,g.outerSources));
z = stagnant_cap_geometry(stagnant_cap_options(struct( ...
    'sourceRequest',50,'phi',0)));
assert(~any(z.cap) && isempty(z.ringOffsets));
zeroBand = stagnant_cap_geometry(stagnant_cap_options(struct( ...
    'sourceRequest',50,'edgeBandAngle',0)));
assert(isempty(zeroBand.ringOffsets) && ~any(zeroBand.ringCollocation));
fprintf('stagnant_cap_geometry self-test passed.\n');
end
