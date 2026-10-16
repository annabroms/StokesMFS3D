function test_ellipsoid_grid_experiment
%TEST_ELLIPSOID_GRID_EXPERIMENT Small shared-grid and Stokes/Laplace regressions.
% Does not invoke either full experiment, sweeps, or figure-generation scripts.

experimentDirectory = fileparts(mfilename('fullpath'));
projectDirectory = fileparts(experimentDirectory);
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
addpath(experimentDirectory,fullfile(projectDirectory,'common','dense'), ...
    fullfile(projectDirectory,'geometry','spheres'), ...
    fullfile(projectDirectory,'geometry','ellipsoids'),fullfile(projectDirectory,'laplace'));

% One tilted prolate pair tests transformations as well as sphere limits.
pair = lap3d_sphere_pair(1,0.05);
angle = 0.37;
rotation = [cos(angle),-sin(angle),0;sin(angle),cos(angle),0;0,0,1];
pair.E1.axes = [1;0.7;0.7];
pair.E2.axes = pair.E1.axes;
tipRadius = 0.7^2;
offset = sqrt(pair.gap*(4*tipRadius+pair.gap))/2;
pair.Pminus = repmat([-offset;0;0],1,721);
pair.Pplus = -pair.Pminus;
pair.E1.Q = rotation;
pair.E2.Q = rotation;
for name = {'p1','p2','n','Pminus','Pplus'}
    pair.(name{1}) = rotation*pair.(name{1});
end
pair.E1.c = rotation*pair.E1.c;
pair.E2.c = rotation*pair.E2.c;
pair.label = 'Tilted tip-to-tip prolate unit test';
opt = lap3d_complete_grid_options(struct('visualise',false,'confirmBeforeSolve',false, ...
    'basicGridDistribution',"sphericaldesign",'basicSourceCount',50, ...
    'basicCollocationCount',200,'enhancingSourceDistribution',"tensor", ...
    'fineSourceResolution',[12,8],'sourceOffset',0.25,'sourceBeta',0.3, ...
    'collocationPatchLayers',8,'testResolution',[12,8], ...
    'testPatchLayers',17,'testPatchAzimuthOversampling',1));
[basic,enhanced] = lap3d_build_two_body_grids(pair,opt);
assert(isequal(basic.testPoints,enhanced.testPoints));
assert(sum(basic.testCounts) == size(basic.testPoints,1));
assert(all(basic.testPatchCounts > 0));
for body = 1:2
    E = pair.(sprintf('E%d',body));
    rows = basic.testBody == body;
    localPoints = (E.Q.'*(basic.testPoints(rows,:).'-E.c))./E.axes;
    assert(max(abs(sum(localPoints.^2,1)-1)) < 1e-12);
    assert(enhanced.sourceCapMetadata{body}.parentNodeCount == 12*8);
end

% Changing enhancing cross-section count changes only the cap, not basic grids.
moreSources = opt;
moreSources.fineSourceResolution = [24,8];
[sameBasic,moreCap] = lap3d_build_two_body_grids(pair,moreSources);
assert(isequal(basic.sources,sameBasic.sources));
assert(isequal(basic.collocation,sameBasic.collocation));
assert(moreCap.sourceCapMetadata{1}.parentNodeCount == 24*8);
assert(size(moreCap.sourceCapBlocks{1},1) > size(enhanced.sourceCapBlocks{1},1));
% Quasiuniform mode also responds to the small-ring azimuthal count floor.
quasi = opt;
quasi.enhancingSourceDistribution = "quasiuniform";
[~,quasiCap] = lap3d_build_two_body_grids(pair,quasi,basic);
quasi.enhancingSourceMinAzimuthNodes = 16;
[~,denserQuasiCap] = lap3d_build_two_body_grids(pair,quasi,basic);
assert(denserQuasiCap.sourceCapMetadata{1}.parentNodeCount > ...
    quasiCap.sourceCapMetadata{1}.parentNodeCount);

% Refining ONLY the check grid must leave source/collocation coordinates unchanged.
finerTests = opt;
finerTests.testResolution = [24,16];
finerTests.testPatchLayers = 25;
finerTests.testPatchAzimuthOversampling = 2;
[fineBasic,fineEnhanced] = lap3d_build_two_body_grids(pair,finerTests);
assert(isequal(basic.sources,fineBasic.sources));
assert(isequal(basic.collocation,fineBasic.collocation));
assert(isequal(enhanced.sources,fineEnhanced.sources));
assert(isequal(enhanced.collocation,fineEnhanced.collocation));
assert(size(fineBasic.testPoints,1) > size(basic.testPoints,1));
testPatch = basic.testGridMetadata.patchMetadata{1};
fitPatch = enhanced.collocationPatchMetadata{1};
assert(testPatch.angularRadius(end) ~= fitPatch.angularRadius(end));

% Batched Laplace evaluation agrees with dense direct evaluation, including inf norms.
laplaceResult = lap3d_run_two_body_grid_experiment(pair,"tip",opt);
laplaceSolution = solve_laplace_capacitance_mfs(enhanced,[1;-1],1e-13,13);
assert(isequal(laplaceResult.enhanced.solution.coefficients,laplaceSolution.coefficients));
densePotential = lap3d_charge_potential_matrix(enhanced.testPoints,enhanced.sources)* ...
    laplaceSolution.coefficients;
assert(norm(densePotential-laplaceSolution.testPotential,inf) < 1e-11);
assert(laplaceSolution.testNearContactMaxResidual == ...
    norm(laplaceSolution.testError(enhanced.testNearContact),inf));
fprintf('Tilted ellipsoid Laplace: independent/near max residual %.3e/%.3e.\n', ...
    laplaceSolution.testMaxResidual,laplaceSolution.testNearContactMaxResidual);

% The shared Stokes driver accepts nonspheres and selected enhancing families.
defaults.enhancingFamilies = 'S';
defaults.translations = [pair.n.';-pair.n.'];
defaults.angularVelocities = zeros(2,3);
defaults.singularityLengthScale = 1;
defaults.evaluationBatchSize = 64;
[controls,remaining] = complete_stokes_two_body_controls( ...
    struct('enhancingFamilies',["D","S"],'testPatchLayers',17),defaults);
assert(strcmp(controls.enhancingFamilies,'SD'));
assert(isequal(remaining,struct('testPatchLayers',17)));
stokesResult = run_stokes_two_body_grid_experiment(pair,"tip",opt,controls);
assert(isequal(stokesResult.enhanced.grid.sources,enhanced.sources));
solution = stokesResult.enhanced.solution;
assert(all(isfinite(solution.coefficients)));
assert(solution.coefficientMax == norm(solution.coefficients,inf));
assert(solution.testMaxResidual == max(solution.testNodeResidual));
assert(solution.testNearContactMaxResidual == max( ...
    solution.testNodeResidual(enhanced.testNearContact)));
assert(solution.unknownCount == 3*(size(basic.sources,1)+ ...
    2*sum(cellfun(@(points) size(points,1),enhanced.sourceCapBlocks))));
fprintf('All small ellipsoid cross-section, evaluation-grid and Stokes-driver tests passed.\n');
end
