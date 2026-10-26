function report = stagnant_cap_validate(solution,validation,opt)
%STAGNANT_CAP_VALIDATE Independent pointwise boundary-condition residuals.
% Each reported RMS is sqrt(mean(pointwise Euclidean residual^2)) over the
% stated sample subset, divided by a fixed physical scale. The sampling
% distribution, rather than surface area, defines this diagnostic RMS.
if nargin == 0
    self_test();
    return
end
opt = stagnant_cap_options(opt);
geometry = solution.grid;
if isfield(solution,'exteriorBasis')
    exteriorBasis = solution.exteriorBasis;
    interiorBasis = solution.interiorBasis;
else
    exteriorBasis = struct('family',"S",'points',geometry.innerSources, ...
        'normalDirections',zeros(0,3),'scale',1);
    interiorBasis = exteriorBasis;
    interiorBasis.points = geometry.outerSources;
end
count = validation.count;
continuity = zeros(count,1);
normal = zeros(count,1);
capOuter = zeros(count,1);
capInner = zeros(count,1);
mobileTraction = zeros(count,1);
surfaceMeridionalVelocity = zeros(count,1);
for first = 1:opt.batchSize:count
    selected = first:min(first+opt.batchSize-1,count);
    x = validation.points(selected,:);
    n = validation.normals(selected,:);
    uo = reshape(stokes_mfs_velocity_matrix(exteriorBasis,x)* ...
        (solution.lambda/opt.muOut),3,[]).';
    ui = reshape(stokes_mfs_velocity_matrix(interiorBasis,x)* ...
        (solution.gamma/opt.muIn),3,[]).';
    to = reshape(stokes_mfs_traction_matrix(exteriorBasis,x,n)* ...
        solution.lambda,3,[]).';
    ti = reshape(stokes_mfs_traction_matrix(interiorBasis,x,n)* ...
        solution.gamma,3,[]).';
    relativeOuter = uo-solution.U(:).';
    relativeInner = ui-solution.U(:).';
    normalOuter = sum(relativeOuter.*n,2);
    normalInner = sum(relativeInner.*n,2);
    jump = to-ti;
    jumpTangent = jump-sum(jump.*n,2).*n;
    continuity(selected) = vecnorm(uo-ui,2,2);
    normal(selected) = abs(normalOuter);
    capOuter(selected) = vecnorm(relativeOuter-normalOuter.*n,2,2);
    % This interior cap trace follows from continuity; it is not a row (32).
    capInner(selected) = vecnorm(relativeInner-normalInner.*n,2,2);
    mobileTraction(selected) = vecnorm(jumpTangent,2,2);
    theta = validation.theta(selected);
    psi = atan2(n(:,2),n(:,1));
    eTheta = [cos(theta).*cos(psi),cos(theta).*sin(psi),-sin(theta)];
    surfaceMeridionalVelocity(selected) = sum(relativeOuter.*eTheta,2);
end
condition = strings(0,1);
region = strings(0,1);
zone = strings(0,1);
nodeCount = zeros(0,1);
minimumEdgeDistance = zeros(0,1);
maximum = zeros(0,1);
rms = zeros(0,1);
names = ["continuity","no penetration","relative cap velocity outer", ...
    "implied interior cap velocity","mobile tangential traction"];
values = {continuity,normal,capOuter,capInner,mobileTraction};
scales = [solution.velocityScale,solution.velocityScale, ...
    solution.velocityScale,solution.velocityScale,solution.tractionScale];
for side = 1:2
    if side==1
        sideMask = validation.cap;
        sideName = "cap";
    else
        sideMask = validation.mobile;
        sideName = "mobile";
    end
    for band = 1:2
        if band==1
            mask = sideMask & validation.nearEdge;
            bandName = "near edge";
        else
            mask = sideMask & ~validation.nearEdge;
            bandName = "elsewhere";
        end
        if ~any(mask), continue; end
        for k = 1:numel(names)
            if (k==3 || k==4) && side~=1 || k==5 && side~=2
                continue
            end
            normalized = values{k}(mask)/scales(k);
            condition(end+1,1) = names(k); %#ok<AGROW>
            region(end+1,1) = sideName; %#ok<AGROW>
            zone(end+1,1) = bandName; %#ok<AGROW>
            nodeCount(end+1,1) = nnz(mask); %#ok<AGROW>
            minimumEdgeDistance(end+1,1) = ...
                min(validation.edgeDistance(mask)); %#ok<AGROW>
            maximum(end+1,1) = max(normalized); %#ok<AGROW>
            rms(end+1,1) = sqrt(mean(normalized.^2)); %#ok<AGROW>
        end
    end
end
report.table = table(condition,region,zone,nodeCount, ...
    minimumEdgeDistance,maximum,rms);
report.point.continuity = continuity;
report.point.normal = normal;
report.point.capOuter = capOuter;
report.point.capInner = capInner;
report.point.mobileTraction = mobileTraction;
report.point.meridionalVelocity = surfaceMeridionalVelocity;
report.minimumEdgeDistance = validation.minimumEdgeDistance;
report.validationCount = count;
end

function self_test()
%SELF_TEST Compare zero residuals and a direct single-Stokeslet trace.
o = stagnant_cap_options(struct('sourceRequest',50,'ringLevels',2, ...
    'validationBaseRequest',250));
g = stagnant_cap_geometry(o);
v = stagnant_cap_validation_grid(o,g.collocationCount,g.ringOffsets(end));
s.grid = g;
s.lambda = zeros(3*g.N,1);
s.gamma = zeros(3*g.N,1);
s.U = zeros(3,1);
s.velocityScale = 1;
s.tractionScale = 1;
r = stagnant_cap_validate(s,v,o);
assert(all(r.table.maximum==0) && all(r.table.rms==0));
s.lambda(1) = 1;
r = stagnant_cap_validate(s,v,o);
x = v.points(1,:);
u = drop_velocity_block(g.innerSources,x,o.muOut)*s.lambda;
assert(abs(r.point.continuity(1)-norm(u))<1e-13);
index = find(v.mobile,1);
n = v.normals(index,:);
t = drop_traction_block(g.innerSources,v.points(index,:),n)*s.lambda;
t = t.'-(t.'*n.')*n;
assert(abs(r.point.mobileTraction(index)-norm(t))<1e-13);
assert(all(isfinite(r.table.maximum)) && max(r.table.maximum)>0);
fprintf('stagnant_cap_validate self-test passed.\n');
end
