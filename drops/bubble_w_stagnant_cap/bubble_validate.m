function report = bubble_validate(solution,validation,opt,cacheFields)
%BUBBLE_VALIDATE Independent exterior boundary-condition residuals.
% Optional cacheFields retains physical fields; existing maxima ordering stays.
if nargin == 0
    bubble_self_test();
    report = [];
    return
end
assert(validation.count==size(validation.points,1));
assert(solution.Uref>0 && solution.tractionScale>0);
if nargin<4, cacheFields = false; end
count = validation.count;
normal = zeros(count,1);
capTangent = zeros(count,1);
mobileTraction = zeros(count,1);
capVelocity = zeros(count,1);
if cacheFields
    report.fields.velocity = zeros(count,3);
    report.fields.relativeVelocity = zeros(count,3);
    report.fields.tangentTraction = zeros(count,3);
end
for first = 1:opt.batchSize:count
    rows = first:min(first+opt.batchSize-1,count);
    x = validation.points(rows,:);
    n = validation.normals(rows,:);
    field = bubble_evaluate(solution,x,opt,n);
    normal(rows) = abs(field.normalVelocity);
    capTangent(rows) = vecnorm(field.tangentVelocity,2,2);
    mobileTraction(rows) = vecnorm(field.tangentTraction,2,2);
    capVelocity(rows) = vecnorm(field.relativeVelocity,2,2);
    if cacheFields
        report.fields.velocity(rows,:) = field.velocity;
        report.fields.relativeVelocity(rows,:) = field.relativeVelocity;
        report.fields.tangentTraction(rows,:) = field.tangentTraction;
    end
end
cap = validation.cap;
mobile = validation.mobile;
maxima = nan(1,4);
if any(cap)
    maxima(1) = max(normal(cap))/solution.Uref;
    maxima(3) = max(capTangent(cap))/solution.Uref;
end
if any(mobile)
    maxima(2) = max(normal(mobile))/solution.Uref;
    maxima(4) = max(mobileTraction(mobile))/solution.tractionScale;
end
report.maxCapVelocity = NaN;
if any(cap), report.maxCapVelocity = max(capVelocity(cap))/solution.Uref; end
report.point.capVelocity = capVelocity;
report.maxima = maxima;
report.conditionNames = ["cap normal velocity", ...
    "mobile normal velocity","cap tangential velocity", ...
    "mobile tangential traction"];
report.point.normal = normal;
report.point.capTangent = capTangent;
report.point.mobileTraction = mobileTraction;
report.validationCount = count;
report.minimumEdgeDistance = validation.minimumEdgeDistance;
report.normalization = [solution.Uref,solution.Uref, ...
    solution.Uref,solution.tractionScale];
end

function bubble_self_test()
% Compare a single Stokeslet trace and verify the four region-specific maxima.
o = bubble_options(struct('sourceRequest',50,'phi',pi/3,'ringLevels',2));
g = bubble_geometry(o);
v = bubble_validation_grid(o,g.collocationCount,g.ringOffsets(end));
s.basis = struct('family',"S",'points',g.sourcePoints, ...
    'normalDirections',zeros(0,3),'scale',1);
s.lambda = zeros(3*g.sourceCount,1);
s.U = zeros(3,1);
s.Uref = 2;
s.tractionScale = 3;
r = bubble_validate(s,v,o);
assert(isequal(r.maxima,zeros(1,4)));
assert(r.maxCapVelocity==0 && ~isfield(r,'fields'));
s.lambda(1) = 1;
r = bubble_validate(s,v,o);
x = v.points(1,:);
n = v.normals(1,:);
u = reshape(stokes_mfs_velocity_matrix(s.basis,x)* ...
    (s.lambda/o.muOut),3,[]).';
assert(abs(r.point.normal(1)-abs(u*n.'))<1e-13);
assert(all(isfinite(r.maxima)) && all(r.maxima>=0));
cached = bubble_validate(s,v,o,true);
assert(isequal(cached.maxima,r.maxima));
assert(norm(cached.point.capVelocity- ...
    sqrt(cached.point.normal.^2+cached.point.capTangent.^2),inf)<1e-13, ...
    'Cap norm mismatch.');
mixed = bubble_options(struct('sourceRequest',50,'phi',pi/3, ...
    'n_level',1,'n_az_source',6,'n_az_colloc',8, ...
    'refinementFamilies',["S","T_t","T_n","T_n2"]));
mixedGrid = bubble_geometry(mixed);
mixedSolution = bubble_solve(mixed,[0;0;-1],mixedGrid);
mixedValidation = bubble_validation_grid( ...
    mixed,mixedGrid.collocationCount,mixedGrid.ringOffsets(end));
mixedReport = bubble_validate(mixedSolution,mixedValidation,mixed);
assert(all(isfinite(mixedReport.maxima)) && ...
    all(mixedReport.maxima>=0) && ...
    mixedReport.validationCount>=5*mixedGrid.collocationCount);
fprintf('bubble_validate self-test passed.\n');
end
