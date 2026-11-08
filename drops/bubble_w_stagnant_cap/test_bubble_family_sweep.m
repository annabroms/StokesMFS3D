function test_bubble_family_sweep()
%TEST_BUBBLE_FAMILY_SWEEP Check paired grids, all kernels, budgets and ranking.
input = struct('sourceRequest',12,'phi',pi/4,'n_level',1,'ringLevels',2, ...
    'n_az_source',6,'n_az_colloc',8,'innerProxyRatio',0.4, ...
    'useColumnScaling',true,'svdTolerance',1e-10,'batchSize',128);
families = {"S","SD","SR","ST_tT_nT_n2","SDRT_tT_nT_n2Q_s","S","S"};
directions = ["none","none","none","none","none","radial","cone_normal"];
outputDirectory = tempname;
mkdir(outputDirectory);
% Remove only this test's newly created temporary output directory on exit.
cleanup = onCleanup(@() rmdir(outputDirectory,'s'));
c = struct('familySets',{families},'derivativeDirections',directions, ...
    'maxMemoryGB',0.05,'maxSvdRows',1500,'maxSvdColumns',200, ...
    'maxCollocationPoints',500,'makePlots',true,'figureVisible','off', ...
    'outputDirectory',outputDirectory);
existingFigures = findall(groot,'Type','figure');
r = bubble_family_sweep(input,c);
newFigures = setdiff(findall(groot,'Type','figure'),existingFigures);
assert(isscalar(newFigures));
assert(numel(findall(newFigures,'Type','axes'))==4);
close(newFigures);
t = r.summary;
assert(all(t.status=="completed") && all(t.forceResidual<1e-7));
assert(all(isfinite(t.validatedMaxima),'all'));
assert(max(abs(t.validatedMaxima-max(t.coarseMaxima,t.fineMaxima)),[],'all')<1e-13);
assert(max(abs(t.overallMaximum-max(t.validatedMaxima,[],2)))<1e-13);
assert(max(abs(t.ratioToStokeslets-t.overallMaximum/t.overallMaximum(1)))<1e-13);
assert(isscalar(unique(t.fitPoints)) && all(t.estimatedGB<=c.maxMemoryGB));
assert(r.fineValidation.count>r.validation.count);
assert(r.fineValidation.minimumEdgeDistance<r.validation.minimumEdgeDistance);
saved = load(fullfile(outputDirectory,'oct29_bubble_families_S_data.mat'));
assert(isequal(saved.solution.grid.points,r.geometries{1}.points));
assert(isequal(saved.check.fineReport.maxima,t.fineMaxima(1,:)));
saved = load(fullfile(outputDirectory,'oct29_bubble_families_summary.mat'));
assert(isequaln(saved.result.summary,t));
assert(isfile(fullfile(outputDirectory,'oct29_bubble_families_summary.csv')));
assert(isfile(fullfile(outputDirectory,'oct29_bubble_families_distance.png')));
for k = 1:numel(families)
    g = r.geometries{k};
    assert(isequal(g.points,r.geometries{1}.points));
    assert(isequal(g.sourcePoints,r.geometries{1}.sourcePoints));
    assert(isequal(g.ringOffsets,r.geometries{1}.ringOffsets));
    assert(all(t.fitPoints(k)>=g.coefficientCount/3));
    assert(isequal(size(r.profiles{k}.maxima),[60 4]));
end
[~,best] = min(t.overallMaximum);
assert(r.lowestSampledCase==t.caseName(best));
if strlength(r.bestResolvedCase)>0
    assert(~t.validationUnresolved(t.caseName==r.bestResolvedCase));
end
% Plan-only mode performs no solves; a column cap excludes enriched bases.
c.planOnly = true;
c.maxSvdColumns = t.columns(1);
plan = bubble_family_sweep(input,c);
assert(plan.summary.status(1)=="planned");
assert(all(plan.summary.status(2:end)=="budget_limit"));
assert(all(isnan(plan.summary.overallMaximum)) && isempty(plan.validation));
c.maxMemoryGB = 1e-8;
c.planOnly = false;
plan = bubble_family_sweep(input,c);
assert(all(plan.summary.status=="budget_limit") && isempty(plan.validation));
saved = load(fullfile(outputDirectory,'oct29_bubble_families_summary.mat'));
assert(all(saved.result.summary.status=="budget_limit"));
fprintf('test_bubble_family_sweep passed: shared grids, all families/derivatives, ranking, plot and budget guards.\n');
end
