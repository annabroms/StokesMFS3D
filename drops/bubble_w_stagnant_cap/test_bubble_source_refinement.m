function test_bubble_source_refinement()
%TEST_BUBBLE_SOURCE_REFINEMENT Check depth-scaled footprints and graph data.
baseOpt = bubble_options(struct('radius',2,'centre',[0.3 -0.2 0.1], ...
    'sourceRequest',12,'phi',pi/3,'n_level',2,'ringLevels',4, ...
    'n_az_source',6,'n_az_colloc',8,'innerSourceStartDistance',0.25, ...
    'sourceClustering',"tapered",'sourceClusterSigma',2, ...
    'collocClusterRatio',0.5));
base = bubble_geometry(baseOpt);
o = baseOpt;
o.useSourceRefinement = true;
o = bubble_options(o);
assert(isequal(o,bubble_options(o)));
g = bubble_geometry(o);
patch = bubble_source_refinement(base,o);
assert(g.sourceRefinementCount==patch.count && patch.count==240);
assert(isequal(g.sourcePoints,base.sourcePoints));
assert(isequal(g.points(1:base.collocationCount,:),base.points));
assert(isequal(g.points(base.collocationCount+1:end,:),patch.points));
assert(all(g.collocationRingIndex(base.collocationCount+1:end)==-1));
assert(g.capCount+g.mobileCount==g.collocationCount);
assert(all(g.sourceRingIndex(patch.sourceIndex)==o.n_level));
radialFactors = patch.edgeDistance./patch.sourceDepth;
assert(all(ismembertol(radialFactors,o.sourceRefinementRadialFactors,1e-12)));
x = patch.points-o.centre;
y = g.sourcePoints(patch.sourceIndex,:)-o.centre;
delta = atan2(x(:,2),x(:,1))-atan2(y(:,2),y(:,1));
delta = atan2(sin(delta),cos(delta));
azimuthFactors = o.radius*sin(patch.theta).*delta./patch.sourceDepth;
assert(all(ismembertol(azimuthFactors,o.sourceRefinementAzimuthFactors,1e-12)));
assert(max(abs(vecnorm(patch.normals,2,2)-1))<1e-13);
scaled = o;
scaled.radius = 2*o.radius;
scaled.centre = 2*o.centre;
scaled.innerSourceStartDistance = 2*o.innerSourceStartDistance;
gScaled = bubble_geometry(bubble_options(scaled));
assert(max(abs(gScaled.points-2*g.points),[],'all')<1e-12);
v = bubble_validation_grid(o,g.collocationCount,g.ringOffsets(end));
fine = bubble_validation_grid(o,2*g.collocationCount,g.ringOffsets(end)/2);
assert(v.count>=5*g.collocationCount && fine.count>v.count);
assert(~any(ismembertol(v.points,g.points,1e-11*o.radius, ...
    'ByRows',true,'DataScale',1)));
audit = bubble_source_refinement(g,o,1);
assert(audit.count>patch.count);
o.useColumnScaling = true;
o.useAreaWeights = true;
s = bubble_solve(o,[0;0;-1],g);
assert(isequal([s.basis.family],["S","S"]));
assert(s.forceResidual<1e-10);
[fitFigure,fitReport] = bubble_plot_collocation_residuals(s,o,'off');
assert(contains(fitFigure.Name,'Collocation points'));
assert(fitReport.validationCount==g.collocationCount);
fitNorm = norm([fitReport.point.normal/s.Uref; ...
    fitReport.point.capTangent(g.cap)/s.Uref; ...
    fitReport.point.mobileTraction(g.mobile)/s.tractionScale]);
assert(abs(fitNorm/sqrt(3*g.collocationCount)-s.fitRMS)<1e-10);
scatters = findall(fitFigure,'Type','Scatter');
pointCount = 0;
for k = 1:numel(scatters)
    pointCount = pointCount+numel(scatters(k).XData);
end
assert(numel(scatters)==4 && pointCount==2*g.collocationCount);
close(fitFigure);
r = bubble_validate(s,v,o);
[f,p] = bubble_plot_residuals_distance(r,v,s,o,'off');
assert(max(abs([p.peakResidual]-r.maxima))<1e-12);
assert(all([p.peakDistanceOverA]>0));
close(f);
bubble_plot_residuals_distance();
fprintf('test_bubble_source_refinement passed: fixed Stokeslets, local scales, independent audit, and plots.\n');
end
