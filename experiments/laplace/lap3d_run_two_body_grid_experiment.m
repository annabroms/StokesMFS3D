function result = lap3d_run_two_body_grid_experiment(pair,caseName,opt)
%LAP3D_RUN_TWO_BODY_GRID_EXPERIMENT Run basic and enhanced capacitance solves.

gridStart = tic;
[basicGrid,enhancedGrid] = lap3d_build_two_body_grids(pair,opt);
gridBuildTime = toc(gridStart);

bodyPotentials = [1;-1];
fprintf('\n%s\n',repmat('=',1,72));
fprintf('Laplace MFS grid experiment: %s\n',pair.label);
fprintf('gap %.12g, prescribed potentials [%g, %g]^T\n', ...
    pair.gap,bodyPotentials);
fprintf('Enhancing Laplace basis: %s; DLP column length scale %g\n', ...
    opt.enhancingLaplaceBasis,opt.laplaceDipoleLengthScale);
fprintf('Column scaling %s; sqrt(surface element) row weighting %s\n', ...
    string(opt.columnScaling),string(opt.surfaceElementPreconditioning));
if opt.columnScaling || opt.surfaceElementPreconditioning
    fprintf('SVD cutoffs/ranks refer to the scaled matrix; residuals remain unweighted.\n');
end
fprintf('%s\n',repmat('=',1,72));
print_focal_diagnostic(basicGrid.focalDiagnostic);
fprintf('\nBasic grid: %s, sources/collocation per body %d/%d\n', ...
    opt.basicGridDistribution,size(basicGrid.sourceBlocks{1},1), ...
    size(basicGrid.collocationBlocks{1},1));
fprintf('\nEnhancing source surface: %s\n',opt.enhancingSourceShape);
if opt.enhancingSourceDistribution == "sphericaldesign"
    fprintf('  requested parent count %d per body, spherical design\n', ...
        opt.enhancingSourceCount);
else
    fprintf('  requested resolution [%d,%d], %s distribution\n', ...
        opt.fineSourceResolution,opt.enhancingSourceDistribution);
end
fprintf('  signed chop offset %+.6g, maximum boundary depth %.6g\n', ...
    opt.sourceChopOffset,opt.sourceOffset+opt.sourceChopOffset);
for body = 1:2
    metadata = enhancedGrid.sourceCapMetadata{body};
    if metadata.shape == "spheroid"
        fprintf(['  body %d: requested beta %.4g, effective beta %.4g, ', ...
            'parent/retained sources %d/%d\n'], ...
            body,opt.sourceBeta,metadata.effectiveBeta, ...
            metadata.parentNodeCount,metadata.retainedCount);
        if metadata.proxyReferenceIsSampled
            fprintf('    nearest sampled proxy-point level %.6g\n', ...
                metadata.nearestProxyLevel);
        else
            fprintf('    continuous proxy radius %.6g; mirrored cap %d\n', ...
                metadata.proxyRadius,metadata.isMirrored);
        end
        fprintf('    width scale %.4g, parent semi-axes [%.6g, %.6g, %.6g]\n', ...
            metadata.widthScale,metadata.parentAxes);
    else
        fprintf(['  body %d: sphere radius %.6g, ', ...
            'parent/retained sources %d/%d\n'], ...
            body,opt.sourceSphereRadius,metadata.parentNodeCount, ...
            metadata.retainedCount);
    end
    fprintf('    retained boundary-depth range [%.6g, %.6g]\n', ...
        metadata.minimumRetainedClearance,metadata.maximumRetainedClearance);
end
fprintf(['\nEnhancing collocation points after duplicate removal: ', ...
    '%d + %d = %d total (excluding basic nodes)\n'], ...
    enhancedGrid.enhancingCollocationCounts, ...
    enhancedGrid.totalEnhancingCollocation);
fprintf('Independent evaluation nodes: %d + %d = %d total (%d near contact)\n', ...
    basicGrid.testCounts,size(basicGrid.testPoints,1),nnz(basicGrid.testNearContact));

% Preview actual dimensions before any dense matrix assembly or SVD.
systemSizes = [size(basicGrid.collocation,1),size(basicGrid.sources,1); ...
    size(enhancedGrid.collocation,1), ...
    size(enhancedGrid.sources,1)+numel(enhancedGrid.dipoleSourceIndices)];
if ~mfs_confirm_system_sizes(systemSizes,["basic Laplace","enhanced Laplace"], ...
        opt.confirmBeforeSolve)
    result = [];
    return
end
fprintf('\nSolving basic grid...\n');
basicSolution = solve_laplace_capacitance_mfs( ...
    basicGrid,bodyPotentials,opt.basicSvdCutoff);
print_solution_summary('basic',basicSolution);

fprintf('\nSolving enhanced grid...\n');
enhancedSolution = solve_laplace_capacitance_mfs( ...
    enhancedGrid,bodyPotentials,opt.enhancedSvdCutoff);
print_solution_summary('enhanced',enhancedSolution);

comparison.relativeChargeDifference = norm( ...
    enhancedSolution.charges-basicSolution.charges)/ ...
    max(norm(enhancedSolution.charges),eps);
comparison.collocationMaxResidualRatio = ...
    enhancedSolution.collocationMaxResidual/ ...
    max(basicSolution.collocationMaxResidual,eps);
comparison.testMaxResidualRatio = enhancedSolution.testMaxResidual/ ...
    max(basicSolution.testMaxResidual,eps);
comparison.testResidualNotMateriallyWorse = ...
    comparison.testMaxResidualRatio <= 1.05;
comparison.nearContactMaxResidualRatio = enhancedSolution.testNearContactMaxResidual/ ...
    max(basicSolution.testNearContactMaxResidual,eps);

fprintf('\nEnhanced/basic relative-residual ratios (inf norm):\n');
fprintf('  collocation %.3e, independent test %.3e\n', ...
    comparison.collocationMaxResidualRatio,comparison.testMaxResidualRatio);
fprintf('  relative change in body charges %.3e\n', ...
    comparison.relativeChargeDifference);

if ~comparison.testResidualNotMateriallyWorse
    warning('lap3d_run_two_body_grid_experiment:GlobalErrorIncrease', ...
        'The enhanced test-node maximum residual exceeds the basic one by over 5%%.');
end

result.caseName = string(caseName);
result.delta = pair.gap;
result.bodyPotentials = bodyPotentials;
result.options = opt;
result.geometry = pair;
result.accumulationPoints = {pair.Pminus,pair.Pplus};
result.gridBuildTime = gridBuildTime;
result.focalDiagnostic = basicGrid.focalDiagnostic;
result.basic.grid = basicGrid;
result.basic.solution = basicSolution;
result.enhanced.grid = enhancedGrid;
result.enhanced.solution = enhancedSolution;
result.enhancingCollocationCounts = enhancedGrid.enhancingCollocationCounts;
result.totalEnhancingCollocation = enhancedGrid.totalEnhancingCollocation;
result.comparison = comparison;
end


function print_focal_diagnostic(diagnostic)
%PRINT_FOCAL_DIAGNOSTIC reports the basic proxy's focal-set relation.

fprintf('\nBasic source proxy and focal set:\n');
fprintf('  focal half-length %.8g, tip clearance %.8g\n', ...
    diagnostic.focalHalfLength,diagnostic.tipClearance);
fprintf('  source offset enclosure margin %.8g: %s\n', ...
    diagnostic.enclosureMargin,string(diagnostic.enclosed));
if diagnostic.proxyIsConfocal
    fprintf('  the spherical inward-offset proxy shares the degenerate focus.\n');
else
    fprintf('  the inward normal-offset proxy is not itself confocal.\n');
end
end


function print_solution_summary(name,sol)
%PRINT_SOLUTION_SUMMARY prints reproducible discretization and error data.

fprintf(['  %s: M=%d, N=%d, cutoff %.1e, retained rank %d, ', ...
    'kept cond %.3e\n'],name,sol.matrixSize(1),sol.matrixSize(2), ...
    sol.svdCutoff,sol.retainedRank,sol.keptConditionNumber);
fprintf('    charges [% .10g, % .10g]^T, coefficient inf norm %.3e\n', ...
    sol.charges,sol.coefficientNorm);
fprintf('    monopole/dipole coefficient maxima %.3e / %.3e (charge units)\n', ...
    sol.monopoleCoefficientMax,sol.dipoleCoefficientMax);
fprintf(['    relative residuals (inf norm): collocation %.3e, ', ...
    'test nodes %.3e, near contact %.3e\n'], ...
    sol.collocationMaxResidual,sol.testMaxResidual,sol.testNearContactMaxResidual);
fprintf('    assembly %.3fs, scaling %.3fs, SVD solve %.3fs\n', ...
    sol.assemblyTime,sol.preconditioningTime,sol.solveTime);
fprintf('    physical evaluation: collocation %.3fs, test %.3fs\n', ...
    sol.collocationEvaluationTime,sol.testEvaluationTime);
end
