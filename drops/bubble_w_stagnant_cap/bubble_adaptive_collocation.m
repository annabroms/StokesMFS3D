function result = bubble_adaptive_collocation(user,controls)
%BUBBLE_ADAPTIVE_COLLOCATION Add large-residual targets under dense SVD limits.
% Keep the Stokeslet source basis fixed. Success requires physical maxima
% below tolerance on fit points, candidates, and two fresh independent audits.
% Finite audits are numerical evidence, not a continuous-surface certificate.
% Refine in complete uniform rings by default, preserving the global seeds.
if nargin<2, controls = struct(); end
defaults = struct('tolerance',1e-2,'maxIterations',15, ...
    'initialCollocationPoints',1200,'maxCollocationPoints',3000, ...
    'maxNewPointsPerIteration',192,'maxSvdColumns',3000, ...
    'selectionMode',"rings",'ringAzimuthCount',64, ...
    'maxSvdRows',9000,'maxMemoryGB',2,'validationBaseRequest',10000, ...
    'edgeAuditDistanceOverA',1e-10,'auditSourceLevels',3, ...
    'useResidualWeights',true,'balanceConditions',true, ...
    'weightDamping',0.5,'maxPointWeight',100, ...
    'maxConditionWeight',30,'stagnationIterations',4, ...
    'minimumImprovement',0.01,'verbose',true);
names = fieldnames(controls);
assert(all(ismember(names,fieldnames(defaults))), ...
    'bubble_adaptive_collocation:UnknownControl','Unknown adaptive control.');
for k = 1:numel(names), defaults.(names{k}) = controls.(names{k}); end
c = defaults;
for name = ["maxIterations","initialCollocationPoints","maxCollocationPoints", ...
        "maxNewPointsPerIteration","maxSvdColumns","maxSvdRows", ...
        "validationBaseRequest","auditSourceLevels","stagnationIterations", ...
        "ringAzimuthCount"]
    validateattributes(c.(name),{'numeric'}, ...
        {'scalar','real','finite','integer','positive'});
end
assert(c.ringAzimuthCount>=4,'Use at least four points per ring.');
c.selectionMode = string(validatestring(c.selectionMode,{'rings','points'}));
if c.selectionMode=="rings"
    assert(c.maxNewPointsPerIteration>=c.ringAzimuthCount, ...
        'An iteration must allow at least one complete refinement ring.');
end
for name = ["tolerance","maxMemoryGB","edgeAuditDistanceOverA", ...
        "maxPointWeight","maxConditionWeight"]
    validateattributes(c.(name),{'numeric'},{'scalar','real','finite','positive'});
end
validateattributes(c.weightDamping,{'numeric'},{'scalar','>',0,'<=',1});
validateattributes(c.minimumImprovement,{'numeric'},{'scalar','>=',0,'<',1});
assert(c.maxPointWeight>=1 && c.maxConditionWeight>=1);
for name = ["useResidualWeights","balanceConditions","verbose"]
    assert(islogical(c.(name)) && isscalar(c.(name)));
end
o = bubble_options(user);
assert(isequal(o.refinementFamilies,"S") && o.derivative_direction=="none", ...
    'bubble_adaptive_collocation:StokesletsOnly','Use Stokeslets only.');
o.maxSvdMemoryGB = min(o.maxSvdMemoryGB,c.maxMemoryGB);
c.maxMemoryGB = o.maxSvdMemoryGB;
full = bubble_geometry(o);
columns = full.coefficientCount;
result = struct('converged',false,'status',"iteration_limit", ...
    'message',"Iteration limit reached before all residuals passed.", ...
    'solution',[],'opt',o,'controls',c,'history',table(), ...
    'fitReport',[],'validation',[],'report',[], ...
    'auditGrids',{{}},'auditReports',{{}});
if columns>c.maxSvdColumns
    result.status = "column_limit";
    result.message = "Source basis exceeds maxSvdColumns; reduce proxy or ring sources.";
    return
end
% Find the largest affordable point count before allocating any dense matrix.
low = 0;
high = min(c.maxCollocationPoints,floor(c.maxSvdRows/3));
while low<high
    trial = ceil((low+high)/2);
    memory = bubble_svd_memory(3*trial,columns);
    if memory.estimatedGB<=c.maxMemoryGB, low = trial; else, high = trial-1; end
end
pointLimit = low;
minimumPoints = max(floor(columns/3)+1,full.baseCollocationCount);
if c.selectionMode=="rings"
    minimumPoints = full.baseCollocationCount+c.ringAzimuthCount* ...
        ceil((minimumPoints-full.baseCollocationCount)/c.ringAzimuthCount);
end
if pointLimit<minimumPoints
    result.status = "memory_limit";
    result.message = "Budget cannot fit the global nodes and an overdetermined SVD.";
    return
end
initial = max(minimumPoints,min([c.initialCollocationPoints,pointLimit, ...
    full.collocationCount]));
base = (1:full.baseCollocationCount).';
extra = (full.baseCollocationCount+1:full.collocationCount).';
if c.selectionMode=="rings"
    angles = uniquetol(full.theta(extra),1e-12,'DataScale',1);
    ringCount = min(numel(angles), ...
        floor((initial-numel(base))/c.ringAzimuthCount));
    if numel(base)+ringCount*c.ringAzimuthCount<minimumPoints
        result.status = "memory_limit";
        result.message = "The seed rings cannot fit an overdetermined solve within the budget.";
        return
    end
    selected = round(linspace(1,numel(angles),ringCount));
    seedPoints = ring_points(angles(selected),o,c.ringAzimuthCount);
    geometry = set_targets(full,[full.points(base,:);seedPoints],o);
else
    selected = round(linspace(1,numel(extra),initial-numel(base)));
    geometry = set_targets(full,full.points([base;extra(selected)],:),o);
end
weights = ones(geometry.collocationCount,3);
baseConditionWeights = o.conditionWeights;
history = zeros(c.maxIterations,14);
bestMaximum = inf;
stagnation = 0;
auditGrids = cell(1,2);
auditReports = cell(1,2);
Fg = (4*pi*o.radius^3/3)*o.densityDifference*o.gravity(:);
for iteration = 1:c.maxIterations
    memory = bubble_svd_memory(3*geometry.collocationCount,columns);
    if c.verbose
        fprintf('\nAdaptive iteration %d: %d points, %d columns, estimate %.3g/%.3g GB.\n', ...
            iteration,geometry.collocationCount,columns,memory.estimatedGB,c.maxMemoryGB);
    end
    solution = bubble_solve(o,Fg,geometry,weights);
    fitGrid = geometry;
    fitGrid.count = geometry.collocationCount;
    fitReport = bubble_validate(solution,fitGrid,o);
    fitErrors = point_errors(fitReport,fitGrid,solution);
    candidate = audit_grid(o,geometry,c,iteration,1);
    candidateReport = bubble_validate(solution,candidate,o);
    fitMax = three_maxima(fitReport);
    candidateMax = three_maxima(candidateReport);
    poolPoints = candidate.points;
    poolErrors = point_errors(candidateReport,candidate,solution);
    % Recheck on fresh denser grids; failed audit points can then become fit points.
    if all([fitMax,candidateMax]<c.tolerance)
        for audit = 1:2
            auditGrids{audit} = audit_grid(o,geometry,c,iteration+0.37*audit,2);
            auditReports{audit} = bubble_validate(solution,auditGrids{audit},o);
            poolPoints = [poolPoints;auditGrids{audit}.points]; %#ok<AGROW>
            poolErrors = [poolErrors; ...
                point_errors(auditReports{audit},auditGrids{audit},solution)]; %#ok<AGROW>
        end
        passed = all([three_maxima(auditReports{1}), ...
            three_maxima(auditReports{2})]<c.tolerance);
        if passed
            result.converged = true;
            result.status = "converged";
            result.message = "All three conditions passed fit, candidate, and two independent audits.";
        end
    end
    worst = max([fitMax,candidateMax,poolErrors(:).']);
    if worst<bestMaximum*(1-c.minimumImprovement)
        bestMaximum = worst;
        stagnation = 0;
    else
        stagnation = stagnation+1;
    end
    history(iteration,:) = [iteration,geometry.collocationCount,columns, ...
        memory.estimatedGB,fitMax,candidateMax,o.conditionWeights,0];
    if c.verbose
        fprintf('Fit maxima: %.3g %.3g %.3g; candidates: %.3g %.3g %.3g.\n', ...
            fitMax,candidateMax);
    end
    if result.converged || iteration==c.maxIterations, break; end
    if stagnation>=c.stagnationIterations && iteration>1 && ...
            history(iteration-1,14)==0
        remaining = pointLimit-geometry.collocationCount;
        if remaining==0 || (c.selectionMode=="rings" && remaining<c.ringAzimuthCount)
            result.status = "point_limit";
            result.message = "Point/memory budget exhausted and reweighting stagnated.";
        else
            result.status = "stagnated";
            result.message = "Residuals stagnated; fixed source basis or SVD truncation may limit accuracy.";
        end
        break
    end
    if c.useResidualWeights
        weights = (weights.*max(fitErrors/c.tolerance,0.1)).^c.weightDamping;
        masks = {true(geometry.collocationCount,1),geometry.cap,geometry.mobile};
        for j = 1:3
            active = masks{j};
            if any(active)
                weights(active,j) = weights(active,j)/mean(weights(active,j));
            end
        end
        weights = min(c.maxPointWeight,max(1/c.maxPointWeight,weights));
    end
    if c.balanceConditions
        conditionError = max(fitMax,candidateMax);
        target = baseConditionWeights.*min(c.maxConditionWeight, ...
            sqrt(max(conditionError/c.tolerance,1)));
        o.conditionWeights = sqrt(o.conditionWeights.*target);
    end
    room = min(c.maxNewPointsPerIteration,pointLimit-geometry.collocationCount);
    additions = select_worst(poolPoints,poolErrors,geometry,o,c,room);
    if ~isempty(additions)
        geometry = set_targets(geometry,[geometry.points;additions],o);
        weights = [weights;ones(size(additions,1),3)]; %#ok<AGROW>
    end
    history(iteration,14) = size(additions,1);
end
% Always return fresh audits, including unsuccessful exits with a solution.
if ~result.converged
    for audit = 1:2
        auditGrids{audit} = audit_grid(o,solution.grid,c,iteration+0.37*audit,2);
        auditReports{audit} = bubble_validate(solution,auditGrids{audit},o);
    end
end
result.solution = solution;
result.opt = o;
result.fitReport = fitReport;
result.validation = auditGrids{2};
result.report = auditReports{2};
result.auditGrids = auditGrids;
result.auditReports = auditReports;
result.history = array2table(history(1:iteration,:),'VariableNames', ...
    {'iteration','points','columns','estimatedGB','fitNormal','fitCapVelocity', ...
    'fitMobileTraction','candidateNormal','candidateCapVelocity', ...
    'candidateMobileTraction','normalWeight','capVelocityWeight', ...
    'mobileTractionWeight','addedPoints'});
if c.verbose
    fprintf('%s: %s\n',result.status,result.message);
    fprintf('Independent audit maxima: %.3g %.3g %.3g.\n',three_maxima(result.report));
end
end

function geometry = set_targets(geometry,points,opt)
% Preserve the source basis and global-node prefix while changing fit targets.
geometry.points = points;
geometry.normals = (points-opt.centre)/opt.radius;
geometry.theta = acos(max(-1,min(1,geometry.normals(:,3))));
geometry.cap = geometry.theta<opt.phi;
geometry.mobile = ~geometry.cap;
geometry.collocationCount = size(points,1);
geometry.capCount = nnz(geometry.cap);
geometry.mobileCount = nnz(geometry.mobile);
geometry.edgeDistance = opt.radius*abs(geometry.theta-opt.phi);
geometry.minimumEdgeDistance = min(geometry.edgeDistance);
geometry.collocationRegion = repmat("mobile",geometry.collocationCount,1);
geometry.collocationRegion(geometry.cap) = "cap";
geometry.collocationRingIndex = -2*ones(geometry.collocationCount,1);
geometry.collocationRingIndex(1:geometry.baseCollocationCount) = 0;
geometry.collocationRingCount = [0 0];
geometry.ringAzimuthCount = 0;
geometry.ringOffsets = [];
geometry.sourceRefinementCount = 0;
end

function grid = audit_grid(opt,geometry,c,phase,refinement)
% Create phase-shifted global/edge samples and source-depth footprint audits.
go = opt.geometryOptions;
go.validationBaseRequest = refinement*c.validationBaseRequest;
if go.phi>0 && go.phi<pi
    go.edgeBandAngle = min([max(go.edgeBandAngle,go.collocAngularExtent), ...
        0.9*go.phi,0.9*(pi-go.phi)]);
end
grid = stagnant_cap_validation_grid(go, ...
    refinement*geometry.collocationCount,2*c.edgeAuditDistanceOverA);
angle = mod(phase*sqrt(2),2*pi);
rotation = [cos(angle),-sin(angle),0;sin(angle),cos(angle),0;0,0,1];
unit = grid.normals*rotation.';
theta = acos(max(-1,min(1,unit(:,3))));
psi = atan2(unit(:,2),unit(:,1));
graded = (grid.baseCount+1:grid.count).';
factor = 0.94+0.04*sin(phase*sqrt(3));
theta(graded) = opt.phi+(theta(graded)-opt.phi)*factor;
unit = [sin(theta).*cos(psi),sin(theta).*sin(psi),cos(theta)];
points = opt.centre+opt.radius*unit;
if opt.n_level>0
    patchOpt = opt;
    patchOpt.sourceRefinementLevels = c.auditSourceLevels;
    auditGeometry = geometry;
    auditGeometry.points = zeros(0,3); % Shift first, then exclude actual fit nodes.
    patch = bubble_source_refinement(auditGeometry,patchOpt,refinement);
    patchTheta = opt.phi+(patch.theta-opt.phi)*(0.96+0.03*cos(phase*sqrt(5)));
    patchPsi = atan2(patch.normals(:,2),patch.normals(:,1));
    patchUnit = [sin(patchTheta).*cos(patchPsi), ...
        sin(patchTheta).*sin(patchPsi),cos(patchTheta)];
    points = [points;opt.centre+opt.radius*patchUnit];
end
% Audit points are never reused as fit points without making a fresh audit.
old = ismembertol(points,geometry.points,1e-12*opt.radius, ...
    'ByRows',true,'DataScale',1);
grid.points = points(~old,:);
grid.normals = (grid.points-opt.centre)/opt.radius;
grid.theta = acos(max(-1,min(1,grid.normals(:,3))));
grid.cap = grid.theta<opt.phi;
grid.mobile = ~grid.cap;
grid.edgeDistance = opt.radius*abs(grid.theta-opt.phi);
grid.minimumEdgeDistance = min(grid.edgeDistance);
grid.nearEdge = grid.edgeDistance<=opt.radius*go.edgeBandAngle;
grid.count = size(grid.points,1);
end

function errors = point_errors(report,grid,solution)
% Express the three physical condition errors in the solver's reference units.
errors = [report.point.normal/solution.Uref, ...
    report.point.capTangent/solution.Uref, ...
    report.point.mobileTraction/solution.tractionScale];
errors(grid.mobile,2) = 0;
errors(grid.cap,3) = 0;
assert(all(isfinite(errors),'all'),'Nonfinite boundary residuals.');
end

function maxima = three_maxima(report)
% Combine cap/mobile normal checks while retaining the other two conditions.
values = report.maxima;
values(isnan(values)) = 0; % A degenerate cap may have an absent region.
maxima = [max(values(1:2)),values(3:4)];
end

function additions = select_worst(points,errors,geometry,opt,c,room)
% Select large errors across conditions, surface-distance bins, and azimuths.
if c.selectionMode=="rings"
    additions = select_rings(points,errors,geometry,opt,c,room);
    return
end
if room<=0
    additions = zeros(0,3);
    return
end
old = ismembertol(points,geometry.points,1e-12*opt.radius, ...
    'ByRows',true,'DataScale',1);
unit = (points-opt.centre)/opt.radius;
theta = acos(max(-1,min(1,unit(:,3))));
cap = theta<opt.phi;
scores = [errors(:,1).*cap,errors(:,1).*(~cap),errors(:,2:3)];
scores(old,:) = 0;
distance = abs(theta-opt.phi);
distanceBin = floor(log(max(distance,1e-14))/log(1.5));
azimuthBin = floor(mod(atan2(unit(:,2),unit(:,1)),2*pi)* ...
    max(16,opt.n_az_source)/(2*pi));
selected = zeros(0,1);
quota = max(1,ceil(room/4));
for condition = 1:4
    [values,order] = sort(scores(:,condition),'descend');
    order = order(values>c.tolerance);
    [~,distinct] = unique([distanceBin(order),azimuthBin(order),cap(order)], ...
        'rows','stable');
    order = order(distinct);
    selected = [selected;order(1:min(quota,numel(order)))]; %#ok<AGROW>
end
selected = unique(selected,'stable');
if numel(selected)<room
    priority = max(scores,[],2);
    priority(selected) = 0;
    [values,order] = sort(priority,'descend');
    order = order(values>c.tolerance);
    selected = [selected;order(1:min(room-numel(selected),numel(order)))];
end
selected = selected(1:min(room,numel(selected)));
additions = points(selected,:);
[~,distinct] = uniquetol(additions,1e-12*opt.radius,'ByRows',true,'DataScale',1);
additions = additions(sort(distinct),:);
end

function additions = select_rings(points,errors,geometry,opt,c,room)
% Use worst azimuthal errors to add complete rings within the point budget.
additions = zeros(0,3);
ringLimit = floor(room/c.ringAzimuthCount);
if ringLimit==0, return; end
unit = (points-opt.centre)/opt.radius;
theta = acos(max(-1,min(1,unit(:,3))));
cap = theta<opt.phi;
scores = errors;
% One representative latitude per logarithmic edge-distance bin and side.
distanceBin = floor(log(max(abs(theta-opt.phi),1e-14))/log(1.5));
[~,~,group] = unique([distanceBin,cap],'rows');
groupCount = max(group);
ringScores = zeros(groupCount,3);
for condition = 1:3
    ringScores(:,condition) = accumarray(group,scores(:,condition), ...
        [groupCount,1],@max,0);
end
% Each representative is the latitude of the worst point in that bin.
[~,order] = sort(max(scores,[],2),'descend');
[~,first] = unique(group(order),'stable');
angles = zeros(groupCount,1);
angles(group(order(first))) = theta(order(first));
existingAngles = geometry.theta(geometry.baseCollocationCount+1:end);
unavailable = ismembertol(angles,existingAngles,1e-12,'DataScale',1) | ...
    sin(angles)<1e-10;
ringScores(unavailable,:) = 0;
for k = 1:ringLimit
    % Cycle conditions; retain the global worst choice when a condition passes.
    condition = mod(k-1,3)+1;
    [value,index] = max(ringScores(:,condition));
    if value<=c.tolerance
        [value,index] = max(max(ringScores,[],2));
    end
    if value<=c.tolerance, break; end
    ring = ring_points(angles(index),opt,c.ringAzimuthCount);
    % Existing nodes can already occupy part of this ring; add only missing nodes.
    old = ismembertol(ring,[geometry.points;additions],1e-12*opt.radius, ...
        'ByRows',true,'DataScale',1);
    additions = [additions;ring(~old,:)]; %#ok<AGROW>
    ringScores(index,:) = 0;
end
end

function points = ring_points(angles,opt,count)
% Generate the same uniformly spaced azimuths for seed and adaptive rings.
psi = 2*pi*(0:count-1).'/count;
points = zeros(numel(angles)*count,3);
for k = 1:numel(angles)
    rows = (k-1)*count+(1:count);
    unit = [sin(angles(k))*cos(psi),sin(angles(k))*sin(psi), ...
        cos(angles(k))*ones(count,1)];
    points(rows,:) = opt.centre+opt.radius*unit;
end
end
