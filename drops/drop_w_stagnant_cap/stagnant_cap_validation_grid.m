function validation = stagnant_cap_validation_grid(opt,maxCollocation,finestOffset)
%STAGNANT_CAP_VALIDATION_GRID Fixed independent design plus finer edge rings.
if nargin == 0
    self_test();
    return
end
opt = stagnant_cap_options(opt);
assert(maxCollocation>0 && maxCollocation==round(maxCollocation));
target = max(opt.validationBaseRequest, ...
    ceil(opt.validationMultiplier*maxCollocation));
try
    [base,~] = stagnant_cap_design(target);
    baseMethod = "spherical design";
catch ME
    if ~strcmp(ME.identifier,'stagnant_cap_design:CatalogueExceeded')
        rethrow(ME);
    end
    % A deterministic quasiuniform sphere preserves independent coverage
    % when the requested validation count exceeds the design catalogue.
    base = fibonacci_sphere(target);
    baseMethod = "Fibonacci sphere";
end
cz = cos(0.37); sz = sin(0.37);
cy = cos(0.23); sy = sin(0.23);
Rz = [cz -sz 0;sz cz 0;0 0 1];
Ry = [cy 0 sy;0 1 0;-sy 0 cy];
base = (Rz*Ry*base.').';
unit = base;
if opt.phi>0 && opt.phi<pi && opt.edgeBandAngle>0
    assert(finestOffset>0 && finestOffset<opt.edgeBandAngle);
    offsets = (finestOffset/2)*2.^(0: ...
        floor(log2(2*opt.edgeBandAngle/finestOffset)));
    offsets = offsets(offsets<=opt.edgeBandAngle);
    nAz = max(24,2*ceil(2*pi*sin(opt.phi)/ ...
        sqrt(4*pi/maxCollocation)));
    rings = zeros(2*numel(offsets)*nAz,3);
    row = 0;
    for k = 1:numel(offsets)
        for side = [-1 1]
            theta = opt.phi+side*offsets(k);
            psi = 2*pi*((0:nAz-1)+0.173+0.31*k)/nAz;
            selected = row+(1:nAz);
            rings(selected,:) = [sin(theta)*cos(psi(:)), ...
                sin(theta)*sin(psi(:)),cos(theta)*ones(nAz,1)];
            row = row+nAz;
        end
    end
    unit = [unit;rings];
end
validation.points = opt.centre+opt.radius*unit;
validation.normals = unit;
validation.theta = acos(max(-1,min(1,unit(:,3))));
validation.cap = validation.theta<opt.phi;
validation.mobile = ~validation.cap;
validation.edgeDistance = opt.radius*abs(validation.theta-opt.phi);
validation.nearEdge = validation.edgeDistance<=opt.radius*opt.edgeBandAngle;
validation.minimumEdgeDistance = min(validation.edgeDistance);
validation.baseCount = size(base,1);
validation.baseMethod = baseMethod;
validation.count = size(unit,1);
assert(validation.count>=opt.validationMultiplier*maxCollocation);
if opt.phi>0 && opt.phi<pi
    assert(validation.minimumEdgeDistance>0);
end
end

function unit = fibonacci_sphere(count)
%FIBONACCI_SPHERE Generate arbitrary-size quasiuniform validation nodes.
k = (0:count-1).';
z = 1-2*(k+0.5)/count;
azimuth = pi*(3-sqrt(5))*k;
radius = sqrt(max(0,1-z.^2));
unit = [radius.*cos(azimuth),radius.*sin(azimuth),z];
end

function self_test()
%SELF_TEST Confirm independent density and closer samples on both sides.
o = stagnant_cap_options(struct('sourceRequest',50));
g = stagnant_cap_geometry(o);
v = stagnant_cap_validation_grid(o,g.collocationCount, ...
    g.ringOffsets(end));
assert(v.count>=5*g.collocationCount);
assert(v.minimumEdgeDistance>0 && v.minimumEdgeDistance<= ...
    o.radius*g.ringOffsets(end)/2+1e-12);
assert(any(v.cap & v.nearEdge) && any(v.mobile & v.nearEdge));
assert(max(abs(vecnorm(v.normals,2,2)-1))<1e-12);
assert(v.baseMethod=="spherical design");
large = o;
largestDesign = size(get_sphdesign(inf),1);
large.validationBaseRequest = largestDesign+1;
f = stagnant_cap_validation_grid(large,g.collocationCount, ...
    g.ringOffsets(end));
assert(f.baseMethod=="Fibonacci sphere" && ...
    f.baseCount==largestDesign+1 && ...
    f.count>=5*g.collocationCount);
assert(max(abs(vecnorm(f.normals,2,2)-1))<1e-12 && ...
    f.minimumEdgeDistance>0 && ...
    any(f.cap & f.nearEdge) && any(f.mobile & f.nearEdge));
fprintf('stagnant_cap_validation_grid self-test passed.\n');
end
