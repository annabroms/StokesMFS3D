function sol = solve_laplace_capacitance_mfs(grid,bodyPotentials,svdCutoff,evaluationBatchSize)
%SOLVE_LAPLACE_CAPACITANCE_MFS Solve a dense exterior Dirichlet MFS problem.
%
%   sol = solve_laplace_capacitance_mfs(grid,bodyPotentials,svdCutoff)
%
% Required grid fields are sources (N-by-3), collocation (M-by-3),
% sourceBody (N-by-1), collocationBody (M-by-1), testPoints (K-by-3), and
% testBody (K-by-1).  The only aggregate errors returned are maximum
% relative residuals at the collocation and independent test nodes:
% ||u-u_prescribed||_inf / ||bodyPotentials||_inf.  coefficientNorm is
% ||coefficients||_inf = max(abs(coefficients)), not the 2-norm.
%
% The overdetermined M-by-N system is solved by a truncated SVD.  A
% singular value is retained when s >= svdCutoff*s(1).  Coefficient order is
% [all source monopoles; enhancing parent-normal dipoles].  Grid dipole
% metadata is optional for legacy monopole-only callers.  DLP columns are
% scaled by dipoleLengthScale L, so physical dipole moments equal L*d*n.
% Each body's net charge sums ONLY monopole coefficients; dipoles have
% zero net charge.  sourceCoefficientMax stores max(|q_j|,|d_j|) per node
% in fitted charge units; coefficientNorm includes BOTH fitted families.
% Optional evaluationBatchSize (default 256) limits memory, not resolution.
% An optional grid.testNearContact mask additionally enables a near-contact
% maximum residual; NaN indicates that no near-contact test nodes were supplied.
%
% Optional grid.columnScaling and grid.surfaceElementPreconditioning default
% to false. Row weights W_ii = sqrt(grid.collocationSurfaceElements(i)) are
% applied to BOTH A and the prescribed data. Column multipliers C_jj are
% 1/||W*A(:,j)||_2. Solve B*z = W*rhs, B = W*A*C, by truncated SVD, then
% return physical coefficients c = C*z. Singular values, retained rank and
% keptConditionNumber refer to B, not the original matrix A. All reported
% residuals, charges and coefficient magnitudes use the unscaled physical
% solution. Area weighting fits an integrated squared residual, not its max.

if nargin == 0
    self_test;
    sol = struct();
    return
end
if nargin < 4
    evaluationBatchSize = 256;
end
validateattributes(evaluationBatchSize,{'numeric'}, ...
    {'real','finite','scalar','integer','positive'});

requiredFields = {'sources','collocation','sourceBody','collocationBody', ...
    'testPoints','testBody'};
for k = 1:numel(requiredFields)
    if ~isfield(grid,requiredFields{k})
        error('solve_laplace_capacitance_mfs:MissingGridField', ...
            'grid.%s is required.',requiredFields{k});
    end
end

% Normalize optional basis metadata once; matrix assembly assumes it exists.
if ~isfield(grid,'dipoleSourceIndices')
    grid.dipoleSourceIndices = zeros(0,1);
    grid.dipoleNormals = zeros(0,3);
    grid.dipoleLengthScale = 1;
end
grid.dipoleSourceIndices = grid.dipoleSourceIndices(:);
booleanNames = {'columnScaling','surfaceElementPreconditioning'};
for k = 1:numel(booleanNames)
    name = booleanNames{k};
    if ~isfield(grid,name)
        grid.(name) = false;
    end
    validateattributes(grid.(name),{'logical','numeric'}, ...
        {'real','scalar','binary'},mfilename,name);
    grid.(name) = logical(grid.(name));
end
sources = grid.sources;
collocation = grid.collocation;
sourceBody = grid.sourceBody(:);
collocationBody = grid.collocationBody(:);
testPoints = grid.testPoints;
testBody = grid.testBody(:);
bodyPotentials = bodyPotentials(:);

N = size(sources,1);
M = size(collocation,1);
K = size(testPoints,1);
P = numel(bodyPotentials);
dipoleCount = numel(grid.dipoleSourceIndices);
unknownCount = N+dipoleCount;
validateattributes(grid.dipoleLengthScale,{'numeric'}, ...
    {'real','finite','scalar','positive'});
assert(isequal(size(grid.dipoleNormals),[dipoleCount,3]), ...
    'Dipole normals must match dipole source indices.');
if dipoleCount > 0
    validateattributes(grid.dipoleSourceIndices,{'numeric'}, ...
        {'real','finite','integer','>=',1,'<=',N});
    assert(numel(unique(grid.dipoleSourceIndices)) == dipoleCount, ...
        'There must be at most one DLP coefficient per physical source node.');
end

if size(sources,2) ~= 3 || size(collocation,2) ~= 3 || size(testPoints,2) ~= 3
    error('solve_laplace_capacitance_mfs:BadCoordinates', ...
        'All coordinate arrays must have three columns.');
end
if numel(sourceBody) ~= N || numel(collocationBody) ~= M || ...
        numel(testBody) ~= K
    error('solve_laplace_capacitance_mfs:BadIndexSize', ...
        'Body-index lengths must match their point arrays.');
end
if M <= unknownCount
    error('solve_laplace_capacitance_mfs:NotOverdetermined', ...
        ['Collocation nodes M=%d must exceed %d unknowns (%d SLP + %d DLP). ', ...
         'Increase basic collocation count/resolution or patch refinement.'], ...
        M,unknownCount,N,dipoleCount);
end
if ~isscalar(svdCutoff) || svdCutoff <= 0 || svdCutoff >= 1
    error('solve_laplace_capacitance_mfs:BadSvdCutoff', ...
        'svdCutoff must lie strictly between zero and one.');
end
if any(sourceBody < 1 | sourceBody > P | sourceBody ~= round(sourceBody)) || ...
        any(collocationBody < 1 | collocationBody > P | ...
            collocationBody ~= round(collocationBody)) || ...
        any(testBody < 1 | testBody > P | testBody ~= round(testBody))
    error('solve_laplace_capacitance_mfs:BadBodyIndex', ...
        'All body indices must be integers between 1 and numel(bodyPotentials).');
end
if any(~isfinite(bodyPotentials))
    error('solve_laplace_capacitance_mfs:BadData', ...
        'Body potentials must be finite.');
end
rowScalingFactors = ones(M,1);
if grid.surfaceElementPreconditioning
    if ~isfield(grid,'collocationSurfaceElements') || ...
            numel(grid.collocationSurfaceElements) ~= M
        error('solve_laplace_capacitance_mfs:MissingSurfaceElements', ...
            'Area preconditioning requires one local surface element per collocation node.');
    end
    validateattributes(grid.collocationSurfaceElements,{'numeric'}, ...
        {'real','finite','vector','positive'});
    rowScalingFactors = sqrt(grid.collocationSurfaceElements(:));
end

assemblyStart = tic;
A = lap3d_mfs_potential_matrix(collocation,grid);
assert(isequal(size(A),[M,unknownCount]),'Laplace basis column count mismatch.');
assemblyTime = toc(assemblyStart);
rhs = bodyPotentials(collocationBody);

preconditioningStart = tic;
scaledRhs = rhs;
if grid.surfaceElementPreconditioning
    A = rowScalingFactors.*A;
    scaledRhs = rowScalingFactors.*rhs;
end
columnScalingFactors = ones(unknownCount,1);
if grid.columnScaling
    columnNorms = vecnorm(A,2,1);
    if any(~isfinite(columnNorms) | columnNorms == 0)
        error('solve_laplace_capacitance_mfs:BadColumnNorm', ...
            'Column scaling requires every weighted column to have positive finite norm.');
    end
    columnScalingFactors = (1./columnNorms).';
    A = A.*columnScalingFactors.';
end
preconditioningTime = toc(preconditioningStart);

solveStart = tic;
[U,S,V] = svd(A,'econ');
singularValues = diag(S);
keep = singularValues >= svdCutoff*singularValues(1);
if ~any(keep)
    error('solve_laplace_capacitance_mfs:EmptyTruncation', ...
        'The SVD cutoff discarded every singular value.');
end
scaledCoefficients = V(:,keep)*((U(:,keep).'*scaledRhs)./singularValues(keep));
coefficients = columnScalingFactors.*scaledCoefficients;
solveTime = toc(solveStart);

% Avoid keeping a second full matrix. Re-evaluate with physical kernels when
% A was scaled, rather than reporting a weighted residual or scaled strengths.
collocationStart = tic;
if grid.columnScaling || grid.surfaceElementPreconditioning
    clear A U S V
    collocationPotential = zeros(M,1);
    for firstNode = 1:evaluationBatchSize:M
        rows = firstNode:min(firstNode+evaluationBatchSize-1,M);
        collocationPotential(rows) = ...
            lap3d_mfs_potential_matrix(collocation(rows,:),grid)*coefficients;
    end
else
    collocationPotential = A*coefficients;
    clear A U S V
end
collocationTime = toc(collocationStart);
collocationError = collocationPotential-rhs;
potentialScale = max(norm(bodyPotentials,inf),eps);
monopoleCoefficients = coefficients(1:N);
dipoleCoefficients = reshape(coefficients(N+1:end),[],1);
charges = accumarray(sourceBody,monopoleCoefficients,[P,1],@sum,0);

testStart = tic;
testPotential = zeros(K,1);
for firstNode = 1:evaluationBatchSize:K
    rows = firstNode:min(firstNode+evaluationBatchSize-1,K);
    testPotential(rows) = ...
        lap3d_mfs_potential_matrix(testPoints(rows,:),grid)*coefficients;
end
testTime = toc(testStart);
testExact = bodyPotentials(testBody);
testError = testPotential-testExact;

sol.coefficients = coefficients;
sol.scaledCoefficients = scaledCoefficients;
sol.columnScaling = grid.columnScaling;
sol.surfaceElementPreconditioning = grid.surfaceElementPreconditioning;
sol.rowScalingFactors = rowScalingFactors;
sol.columnScalingFactors = columnScalingFactors;
sol.svdMatrixDescription = "B = diag(rowScalingFactors)*A*diag(columnScalingFactors)";
sol.monopoleCoefficients = monopoleCoefficients;
sol.dipoleCoefficients = dipoleCoefficients;
sol.dipoleSourceIndices = grid.dipoleSourceIndices;
sol.dipoleLengthScale = grid.dipoleLengthScale;
sol.physicalDipoleMoments = zeros(dipoleCount,3);
sol.coefficientSourceIndices = [(1:N).';grid.dipoleSourceIndices];
sol.coefficientBody = sourceBody(sol.coefficientSourceIndices);
sol.coefficientFamily = [repmat("SLP",N,1);repmat("DLP",dipoleCount,1)];
sol.sourceCoefficientMax = abs(monopoleCoefficients);
if dipoleCount > 0
    sol.physicalDipoleMoments = grid.dipoleLengthScale* ...
        dipoleCoefficients.*grid.dipoleNormals;
    sol.sourceCoefficientMax(grid.dipoleSourceIndices) = max( ...
        sol.sourceCoefficientMax(grid.dipoleSourceIndices),abs(dipoleCoefficients));
end
sol.monopoleCoefficientMax = norm(monopoleCoefficients,inf);
sol.dipoleCoefficientMax = norm(dipoleCoefficients,inf);
sol.charges = charges;
sol.svdCutoff = svdCutoff;
sol.singularValues = singularValues;
sol.retainedRank = nnz(keep);
sol.keptConditionNumber = singularValues(1)/singularValues(find(keep,1,'last'));
sol.collocationPotential = collocationPotential;
sol.collocationExact = rhs;
sol.collocationError = collocationError;
sol.collocationMaxResidual = norm(collocationError,inf)/potentialScale;
sol.testPotential = testPotential;
sol.testError = testError;
sol.testMaxResidual = norm(testError,inf)/potentialScale;
sol.testNearContactMaxResidual = NaN;
if isfield(grid,'testNearContact') && any(grid.testNearContact)
    sol.testNearContactMaxResidual = ...
        norm(testError(grid.testNearContact),inf)/potentialScale;
end
sol.coefficientNorm = norm(coefficients,inf);
sol.assemblyTime = assemblyTime;
sol.preconditioningTime = preconditioningTime;
sol.solveTime = solveTime;
sol.collocationEvaluationTime = collocationTime;
sol.testEvaluationTime = testTime;
sol.sourceCount = N;
sol.unknownCount = unknownCount;
sol.matrixSize = [M,unknownCount];
end


function self_test()
%SELF_TEST checks the kernel normalization through a unit-sphere capacity.

radius = 1.7;
collocation = radius*[1,0,0; -1,0,0; 0,1,0; 0,-1,0; 0,0,1; 0,0,-1];
testPoints = radius*[1,1,1; -1,1,1; 1,-1,1; 1,1,-1];
testPoints = radius*testPoints./vecnorm(testPoints,2,2);

grid.sources = [0,0,0];
grid.collocation = collocation;
grid.sourceBody = 1;
grid.collocationBody = ones(size(collocation,1),1);
grid.testPoints = testPoints;
grid.testBody = ones(size(testPoints,1),1);

sol = solve_laplace_capacitance_mfs(grid,1,1e-14);
assert(abs(sol.charges-4*pi*radius)/(4*pi*radius) < 100*eps);
assert(sol.testMaxResidual < 100*eps);

% Two nonzero coefficients distinguish the max norm from the 2-norm.
shift = [10,0,0];
grid.sources = [0,0,0;shift];
grid.sourceBody = [1;2];
grid.collocation = [collocation;collocation+shift];
grid.collocationBody = [ones(size(collocation,1),1); ...
    2*ones(size(collocation,1),1)];
grid.testPoints = [testPoints;testPoints+shift];
grid.testBody = [ones(size(testPoints,1),1); ...
    2*ones(size(testPoints,1),1)];
bodyPotentials = [2;-1];
sol = solve_laplace_capacitance_mfs(grid,bodyPotentials,1e-14);
assert(sol.coefficientNorm == max(abs(sol.coefficients)));
assert(sol.coefficientNorm < norm(sol.coefficients,2));
assert(sol.collocationMaxResidual == ...
    max(abs(sol.collocationError))/max(abs(bodyPotentials)));
assert(sol.testMaxResidual == ...
    max(abs(sol.testError))/max(abs(bodyPotentials)));
fprintf('solve_laplace_capacitance_mfs self-test passed.\n');
end
