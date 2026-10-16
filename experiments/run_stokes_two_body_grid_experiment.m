function result = run_stokes_two_body_grid_experiment(pair,caseName,opt,controls)
%RUN_STOKES_TWO_BODY_GRID_EXPERIMENT Shared basic/enhanced resistance driver.
% Geometry grids are identical to Laplace; solve prescribed rigid velocities
% with free-space viscosity-one kernels and separate relative SVD cutoffs.

gridStart = tic;
[basicGrid,enhancedGrid] = lap3d_build_two_body_grids(pair,opt);
gridBuildTime = toc(gridStart);
fprintf('\nStokes resistance: %s, gap=%g, viscosity=1\n',pair.label,pair.gap);
fprintf('Enhancing families: %s; basic sources: S only\n',controls.enhancingFamilies);
fprintf('Prescribed translations (body rows):\n');
disp(controls.translations);
fprintf('Basic sources/collocation per body: %d/%d\n', ...
    size(basicGrid.sourceBlocks{1},1),size(basicGrid.collocationBlocks{1},1));
for body = 1:2
    metadata = enhancedGrid.sourceCapMetadata{body};
    fprintf('Body %d: parent/retained sources %d/%d; requested/effective beta %g/%g\n', ...
        body,metadata.parentNodeCount,metadata.retainedCount, ...
        opt.sourceBeta,metadata.effectiveBeta);
    fprintf('  parent semi-axes [%g,%g,%g]; width %g; chop offset %+.6g\n', ...
        metadata.parentAxes,opt.sourceWidthScale,opt.sourceChopOffset);
end
fprintf('Total enhancing collocation: %d + %d = %d (excluding basic)\n', ...
    enhancedGrid.enhancingCollocationCounts,enhancedGrid.totalEnhancingCollocation);
fprintf('Independent evaluation nodes: %d + %d = %d total (%d near contact)\n', ...
    basicGrid.testCounts,size(basicGrid.testPoints,1),nnz(basicGrid.testNearContact));
fprintf('Enhanced velocity equations/unknowns: %d/%d\n', ...
    3*size(enhancedGrid.collocation,1), ...
    3*(size(basicGrid.sources,1)+numel(controls.enhancingFamilies)* ...
    sum(cellfun(@(points) size(points,1),enhancedGrid.sourceCapBlocks))));

centres = [pair.E1.c.';pair.E2.c.'];
% Preserve the operator; obtain consent before any dense kernel allocation.
capCount = sum(cellfun(@(points) size(points,1),enhancedGrid.sourceCapBlocks));
systemSizes = [3*size(basicGrid.collocation,1),3*size(basicGrid.sources,1); ...
    3*size(enhancedGrid.collocation,1), ...
    3*(size(basicGrid.sources,1)+numel(controls.enhancingFamilies)*capCount)];
if ~mfs_confirm_system_sizes(systemSizes,["basic Stokes","enhanced Stokes"], ...
        opt.confirmBeforeSolve)
    result = [];
    return
end
basicSolution = solve_stokes_two_body_mfs(basicGrid,centres,controls,opt.basicSvdCutoff);
print_solution('basic',basicSolution);
enhancedSolution = solve_stokes_two_body_mfs(enhancedGrid,centres,controls,opt.enhancedSvdCutoff);
print_solution('enhanced',enhancedSolution);
comparison.collocationMaxResidualRatio = enhancedSolution.collocationMaxResidual/ ...
    max(basicSolution.collocationMaxResidual,eps);
comparison.testMaxResidualRatio = enhancedSolution.testMaxResidual/ ...
    max(basicSolution.testMaxResidual,eps);
comparison.nearContactMaxResidualRatio = enhancedSolution.testNearContactMaxResidual/ ...
    max(basicSolution.testNearContactMaxResidual,eps);
comparison.relativeForceDifference = max(abs( ...
    enhancedSolution.forcesOnFluid-basicSolution.forcesOnFluid),[],'all')/ ...
    max(max(abs(enhancedSolution.forcesOnFluid),[],'all'),eps);
fprintf('\nEnhanced/basic max-residual ratios: collocation %.3e, test %.3e, near %.3e\n', ...
    comparison.collocationMaxResidualRatio,comparison.testMaxResidualRatio, ...
    comparison.nearContactMaxResidualRatio);
fprintf('Relative force change (max norm): %.3e\n',comparison.relativeForceDifference);

result.caseName = string(caseName);
result.delta = pair.gap;
result.options = opt;
result.stokesControls = controls;
result.geometry = pair;
result.accumulationPoints = {pair.Pminus,pair.Pplus};
result.gridBuildTime = gridBuildTime;
result.basic.grid = basicGrid;
result.basic.solution = basicSolution;
result.enhanced.grid = enhancedGrid;
result.enhanced.solution = enhancedSolution;
result.enhancingCollocationCounts = enhancedGrid.enhancingCollocationCounts;
result.totalEnhancingCollocation = enhancedGrid.totalEnhancingCollocation;
result.comparison = comparison;
end


function print_solution(name,solution)
% Report both solves identically, including force signs and max-norm errors.
fprintf('\n%s: relative velocity residuals in max norm\n',name);
fprintf('  collocation %.3e; independent %.3e; near contact %.3e\n', ...
    solution.collocationMaxResidual,solution.testMaxResidual, ...
    solution.testNearContactMaxResidual);
fprintf('  max coefficient %.3e; retained rank %d/%d; relative SVD cutoff %.1e\n', ...
    solution.coefficientMax,solution.retainedRank,solution.unknownCount,solution.svdCutoff);
fprintf('  assembly / SVD solve / evaluation: %.3g / %.3g / %.3g seconds\n', ...
    solution.assemblyTime,solution.solveTime,solution.evaluationTime);
fprintf('  force on fluid (body rows); force on particles is its negative:\n');
disp(solution.forcesOnFluid);
fprintf('  torque on fluid about each body centre:\n');
disp(solution.torquesOnFluid);
disp(solution.familySummary);
end
