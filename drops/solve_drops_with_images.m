function [U,lambda,gamma,grid,info] = solve_drops_with_images(q,a,S1,e,opt)
%SOLVE_DROPS_WITH_IMAGES Active spherical drops with mixed image sources.
% [U,lambda,gamma,grid,info] = solve_drops_with_images(q,a,S1,e,opt)
% Rows are body-major [continuity;normal velocity;tangential traction].
% Columns per body are [exterior-flow coefficients;interior-flow coefficients].
if nargin == 0
    self_test(); U=[]; lambda={}; gamma={}; grid=struct(); info=struct(); return
end
P = size(q,1);
validateattributes(q,{'numeric'},{'2d','ncols',3,'finite','real'});
if isscalar(a), a = repmat(a,P,1); end
a = a(:);
validateattributes(a,{'numeric'},{'numel',P,'positive','finite'});
if isscalar(S1), S1 = repmat(S1,P,1); end
S1 = S1(:);
validateattributes(S1,{'numeric'},{'numel',P,'finite','real'});
if isequal(size(e),[1 3]), e = repmat(e,P,1); end
if ~isequal(size(e),[P 3]) || any(vecnorm(e,2,2)==0)
    error('solve_drops_with_images:BadDirections','Supply one nonzero direction per drop.');
end
e = e./vecnorm(e,2,2);
if nargin < 5, opt = struct(); end
opt = drop_image_options(opt,P,a);

% Construct and, if requested, display geometry before any factorisation.
setupClock = tic;
grid = drop_image_discretization(q,a,opt);
if opt.plot_grid
    drop_image_plot_grid(grid);
end
M = arrayfun(@(b) size(b.collocation,1),grid.body);
factors = cell(P,1);
newFactorCount = 0;
factorizationTime = 0;
if opt.solverMode=="denseSVD"
    if P~=2
        error('solve_drops_with_images:DenseTwoOnly', ...
            'Dense SVD is provided for exactly two drops.');
    end
    % Dense assembly needs the local blocks but no local pseudoinverses.
    for i = 1:P
        factors{i} = drop_image_scaled_self_block(grid.body(i), ...
            opt.interiorViscosity(i),opt.completionScale(i),opt);
    end
else
    keys = cell(P,1);
    for i = 1:P
        body = grid.body(i);
        key = body;
        % Round only the translation-invariant cache key; nodes stay exact.
        key.collocation = round(body.collocation-q(i,:),12);
        for k = 1:numel(key.basis)
            key.basis(k).points = round(key.basis(k).points-q(i,:),12);
            key.basis(k).normalDirections = ...
                round(key.basis(k).normalDirections,12);
        end
        for k = 1:numel(key.interiorBasis)
            key.interiorBasis(k).points = ...
                round(key.interiorBasis(k).points-q(i,:),12);
            key.interiorBasis(k).normalDirections = ...
                round(key.interiorBasis(k).normalDirections,12);
        end
        key = {key.basis,key.interiorBasis,key.collocation, ...
            round(key.normals,12),round(key.rowScale,12), ...
            opt.interiorViscosity(i),opt.completionScale(i), ...
            opt.velocityScale,opt.tractionScale,opt.tractionRowWeight, ...
            opt.scale_columns,opt.svdCutoff};
        reused = false;
        for j = 1:i-1
            if isequal(key,keys{j})
                factors{i} = factors{j};
                reused = true;
                break
            end
        end
        if ~reused
            factorClock = tic;
            factors{i} = drop_image_one_body_preconditioner(body, ...
                opt.interiorViscosity(i),opt.completionScale(i),opt);
            factorizationTime = factorizationTime+toc(factorClock);
            newFactorCount = newFactorCount+1;
        end
        keys{i} = key;
    end
end

% The prescribed active stress appears only in tangent-traction rows.
rhs = zeros(6*sum(M),1);
row0 = 0;
for i = 1:P
    n = grid.body(i).normals;
    d = n*e(i,:).';
    t = S1(i)*[d.*(factors{i}.tangentOne*e(i,:).'), ...
        d.*(factors{i}.tangentTwo*e(i,:).')];
    t = (grid.body(i).rowScale.*t)* ...
        opt.tractionRowWeight/opt.tractionScale;
    rhs(row0+4*M(i)+(1:2*M(i))) = reshape(t.',[],1);
    row0 = row0+6*M(i);
end
setupTime = toc(setupClock);
matvecTime = NaN;
if opt.solverMode=="oneBody"
    matvecClock = tic;
    drop_image_matvec(zeros(size(rhs)),grid,factors,opt);
    matvecTime = toc(matvecClock);
end

% Solve either the identity-diagonal boundary-data system or a dense reference.
solveClock = tic;
if opt.solverMode=="oneBody"
    if norm(rhs)==0
        eta = zeros(size(rhs)); iters = 0; history = 0;
    else
        [eta,iters,history] = helsing_gmres( ...
            @(v) drop_image_matvec(v,grid,factors,opt),rhs,numel(rhs), ...
            opt.maxIterations,opt.gmresTolerance,0,0);
    end
    z = zeros(sum(cellfun(@(f) size(f.selfMatrix,2),factors)),1);
    row0 = 0; col0 = 0;
    discarded = zeros(size(rhs));
    for i = 1:P
        rows = row0+(1:6*M(i));
        cols = col0+(1:size(factors{i}.selfMatrix,2));
        z(cols) = applyOneBodyPinv(eta(rows),factors{i});
        L = factors{i}.leftTranspose;
        discarded(rows) = eta(rows)-L.'*(L*eta(rows));
        row0 = row0+6*M(i); col0 = col0+numel(cols);
    end
    regularizedResidual = norm(drop_image_matvec(eta,grid,factors,opt)-rhs)/ ...
        max(norm(rhs),realmin);
    rank = cellfun(@(f) f.rank,factors);
else
    cols = cellfun(@(f) size(f.selfMatrix,2),factors);
    A = zeros(6*sum(M),sum(cols));
    r0 = 0;
    for i = 1:P
        c0 = 0;
        for j = 1:P
            rows = r0+(1:6*M(i)); columns = c0+(1:cols(j));
            if i==j
                A(rows,columns) = factors{i}.selfMatrix;
            else
                block = drop_image_interaction_block( ...
                    grid.body(j),grid.body(i),factors{i}.tangentOne, ...
                    factors{i}.tangentTwo,opt.interiorViscosity(i));
                A(rows,columns) = factors{i}.rowScale.*block;
            end
            c0 = c0+cols(j);
        end
        r0 = r0+6*M(i);
    end
    [m,n] = size(A);
    storage = whos('A');
    r = min(m,n);
    factorMiB = (m+n+r)*r*(storage.bytes/numel(A))/2^20;
    fprintf(['drop global SVD: %d x %d, matrix %.2f MiB, ', ...
        'economy factors %.2f MiB estimated\n'], ...
        m,n,storage.bytes/2^20,factorMiB);
    [L,D,R] = svd(A,'econ');
    s = diag(D);
    keep = s>opt.svdCutoff*s(1);
    rank = nnz(keep);
    fprintf('drop global retained rank %d/%d at relative cutoff %.3e\n', ...
        rank,numel(s),opt.svdCutoff);
    z = R(:,keep)*((L(:,keep)'*rhs)./s(keep));
    iters = 0; history = [];
    discarded = rhs-L(:,keep)*(L(:,keep)'*rhs);
    regularizedResidual = NaN;
end
solveTime = toc(solveClock);

% Recover physical strengths and translations from the common coefficient layout.
U = zeros(P,3); lambda = cell(P,1); gamma = cell(P,1);
forceSum = zeros(P,1);
col0 = 0;
for i = 1:P
    H = factors{i}.H;
    ne = size(H,2);
    ni = size(factors{i}.selfMatrix,2)-ne;
    alpha = z(col0+(1:ne));
    gamma{i} = z(col0+ne+(1:ni));
    lambda{i} = alpha-factors{i}.C*(H*alpha);
    U(i,:) = (-opt.completionScale(i)*H*alpha).';
    forceSum(i) = norm(H*lambda{i});
    col0 = col0+ne+ni;
end
originalVector = drop_image_apply_A(z,grid,factors,opt)-rhs;
originalResidual = norm(originalVector)/max(norm(rhs),realmin);
info.solverMode = opt.solverMode;
info.retainedRanks = rank;
info.svdCutoff = opt.svdCutoff;
info.iterations = iters;
info.residualHistory = history;
info.regularizedResidual = regularizedResidual;
info.originalResidual = originalResidual;
info.converged = NaN;
if opt.solverMode=="oneBody"
    info.converged = regularizedResidual<=opt.gmresTolerance;
end
info.residualIdentityError = NaN;
if opt.solverMode=="oneBody"
    regularizedVector = drop_image_matvec(eta,grid,factors,opt)-rhs;
    info.residualIdentityError = norm( ...
        originalVector-(regularizedVector-discarded))/ ...
        max(1,norm(originalVector));
end
info.discardedRangeNorm = norm(discarded)/max(norm(rhs),realmin);
info.forceSumNorm = forceSum;
info.newFactorCount = newFactorCount;
info.factorizationTime = factorizationTime;
info.setupTime = setupTime;
info.matvecTime = matvecTime;
info.solveTime = solveTime;
info.systemSize = [6*sum(M),numel(z)];
info.rhsScaled = rhs;
info.options = opt;
end

function self_test()
% Check that disabling images reproduces the existing basic drop solver.
root = fileparts(fileparts(mfilename('fullpath')));
oldPath = path; cleanup = onCleanup(@() path(oldPath));
addpath(genpath(root));
q = [-1.6 0 0;1.6 0 0]; a = ones(2,1);
opt = struct('N',30,'collocationRatio',2,'useImages',false, ...
    'svdCutoff',1e-12,'gmresTolerance',1e-9,'maxIterations',50);
[U,lambda,gamma,g,info] = solve_drops_with_images(q,a,1,[0 0 1],opt);
old = struct('interiorViscosity',2,'completionScale',1/(8*pi), ...
    'velocityScale',1,'tractionScale',1,'tractionRowWeight',1, ...
    'oneBodySvdTolerance',1e-12,'krylovTolerance',1e-9, ...
    'maxIterations',50,'scale_columns',false,'interactionBatchSize',128);
inner = cell(2,1); outer = cell(2,1); surface = cell(2,1); normals = cell(2,1);
for i = 1:2
    inner{i} = g.body(i).basis(1).points;
    outer{i} = g.body(i).interiorBasis(1).points;
    surface{i} = g.body(i).collocation;
    normals{i} = g.body(i).normals;
end
[U0,~,gamma0,lambda0,oldInfo] = solve_drops_preconditioned( ...
    inner,outer,surface,normals,1,repmat([0 0 1],2,1),old);
assert(norm(U-U0,inf)<1e-8 && ...
    abs(info.originalResidual-oldInfo.originalResidual)<1e-8);
assert(info.residualIdentityError<1e-10);
for i = 1:2
    assert(norm(lambda{i}-lambda0{i},inf)<1e-7);
    assert(norm(gamma{i}-gamma0{i},inf)<1e-7);
end
near = [-1.025 0 0;1.025 0 0];
opt.N_image = 1;
opt.imageFamilies = "SRTD";
opt.mirror_ext = true;
opt.useImages = true;
opt.capPointCounts = [12 12];
[~,~,~,~,iterative] = solve_drops_with_images( ...
    near,a,1,[0 0 1],opt);
opt.solverMode = "denseSVD";
[~,~,~,~,dense] = solve_drops_with_images(near,a,1,[0 0 1],opt);
assert(isequal(iterative.systemSize,dense.systemSize));
assert(dense.originalResidual<=iterative.originalResidual+1e-9);
assert(dense.newFactorCount==0 && dense.factorizationTime==0);
fprintf('solve_drops_with_images self-test passed.\n');
end
