function test_bubble_adaptive_collocation()
%TEST_BUBBLE_ADAPTIVE_COLLOCATION Verify adaptation, weights, and size guards.
input = struct('sourceRequest',12,'phi',pi/4,'n_level',1, ...
    'ringLevels',2,'n_az_source',6,'n_az_colloc',8, ...
    'innerProxyRatio',0.4,'useColumnScaling',true, ...
    'useConditionCountScaling',true);
c = struct('tolerance',1e-12,'maxIterations',2, ...
    'initialCollocationPoints',40,'maxCollocationPoints',60, ...
    'maxNewPointsPerIteration',8,'maxSvdRows',180,'maxSvdColumns',100, ...
    'selectionMode',"rings",'ringAzimuthCount',8, ...
    'maxMemoryGB',0.02,'validationBaseRequest',100,'auditSourceLevels',1, ...
    'verbose',false);
r = bubble_adaptive_collocation(input,c);
assert(~r.converged && height(r.history)==2);
assert(r.history.addedPoints(1)>0 && r.history.points(2)>r.history.points(1));
assert(all(r.history.points<=c.maxCollocationPoints));
assert(all(3*r.history.points<=c.maxSvdRows));
assert(all(r.history.estimatedGB<=c.maxMemoryGB));
assert(r.solution.forceResidual<1e-9);
g = r.solution.grid;
assert(isequal([r.solution.basis.family],["S","S"]));
% Both seed and adaptive targets beyond the global design form complete rings.
ringNodes = g.normals(g.baseCollocationCount+1:end,:);
assert(mod(size(ringNodes,1),c.ringAzimuthCount)==0);
psi = 2*pi*(0:c.ringAzimuthCount-1).'/c.ringAzimuthCount;
for first = 1:c.ringAzimuthCount:size(ringNodes,1)
    ring = ringNodes(first:first+c.ringAzimuthCount-1,:);
    assert(max(abs(ring(:,3)-ring(1,3)))<1e-12);
    radial = sqrt(1-ring(1,3)^2);
    assert(max(abs(ring(:,1:2)-radial*[cos(psi),sin(psi)]),[],'all')<1e-11);
end
assert(mod(r.history.addedPoints(1),c.ringAzimuthCount)==0);
for audit = 1:2
    assert(~any(ismembertol(r.auditGrids{audit}.points,g.points, ...
        1e-12,'ByRows',true,'DataScale',1)));
end
% A remainder smaller than a full ring must remain unused.
limited = c;
limited.maxCollocationPoints = r.history.points(1)+c.ringAzimuthCount-1;
blocked = bubble_adaptive_collocation(input,limited);
assert(all(blocked.history.addedPoints==0));
assert(all(blocked.history.points==r.history.points(1)));
% Individual-point mode is still available for comparisons.
individual = c;
individual.selectionMode = "points";
pointResult = bubble_adaptive_collocation(input,individual);
assert(pointResult.history.addedPoints(1)>0 && ~pointResult.converged);
% A tiny budget must refuse BEFORE attempting matrix assembly or an SVD.
limited = c;
limited.maxMemoryGB = 1e-8;
blocked = bubble_adaptive_collocation(input,limited);
assert(blocked.status=="memory_limit" && isempty(blocked.solution));
limited = c;
limited.maxSvdColumns = 10;
blocked = bubble_adaptive_collocation(input,limited);
assert(blocked.status=="column_limit" && isempty(blocked.solution));
o = bubble_options(input);
o.maxSvdMemoryGB = 1e-8;
try
    bubble_solve(o,[0;0;-1]);
    error('test_bubble_adaptive_collocation:MissingGuard','Budget was ignored.');
catch ME
    assert(strcmp(ME.identifier,'bubble_solve:MemoryBudget'));
end
% Independently reconstruct row weights in the three projected row blocks.
o = bubble_options(input);
o.conditionWeights = [2 3 4];
g = bubble_geometry(o);
weights = ones(g.collocationCount,3);
weights(:,1) = linspace(0.2,2,g.collocationCount).';
s = bubble_solve(o,[0;0;-1],g,weights);
expected = [2*sqrt(weights(:,1)/g.collocationCount)/s.Uref; ...
    repelem(3*sqrt(weights(g.cap,2)/g.capCount)/s.Uref,2); ...
    repelem(4*sqrt(weights(g.mobile,3)/g.mobileCount)/s.tractionScale,2)];
assert(norm(s.rowScale-expected,inf)<1e-12);
% A smooth endpoint case must genuinely pass the default 1e-2 target.
smooth = struct('sourceRequest',50,'phi',0,'n_level',0,'ringLevels',0, ...
    'innerProxyRatio',0.2,'useColumnScaling',true, ...
    'useConditionCountScaling',true);
passed = c;
passed.tolerance = 1e-2;
passed.initialCollocationPoints = 90;
passed.maxCollocationPoints = 120;
passed.maxSvdRows = 360;
passed.maxSvdColumns = 200;
successful = bubble_adaptive_collocation(smooth,passed);
assert(successful.converged && successful.status=="converged");
for audit = 1:2
    maxima = successful.auditReports{audit}.maxima;
    assert(all(maxima(isfinite(maxima))<passed.tolerance));
end
fprintf('test_bubble_adaptive_collocation passed: uniform rings, point mode, SVD limits, weighting, independent success/failure.\n');
end
