function [U,iters,lambda_norm,uerr,details] = ...
    solve_mobility_images(q,Fvec,opt)
%SOLVE_MOBILITY_IMAGES Solve unit-sphere mobility by recompleted image MFS.
%   [U,iters,lambda_norm,uerr,details] = solve_mobility_images(q,Fvec,opt)
%   solve_mobility_images() runs a small resistance-to-mobility check.
%
% q is P-by-3; Fvec and U are 6P-by-1, with each body ordered
% [force;torque] and [translation;angular velocity], respectively.
% Fvec uses the source-strength, force-on-fluid sign of the image resistance
% solver. opt shares its geometry, row weights, SVD cutoff, and GMRES
% controls; useRowScaling=false selects identity row scaling, while
% solverMode is "oneBody" (default) or "denseSVD". details.grid and
% details.lambda contain the complete basis and physical strengths for
% evaluating the flow with direct_image_velocity or the dense kernel.
%
% The unprojected unknown z stores rigid motion through U_k=-H_k*z_k.
% Its projected part has zero total force and torque across all families;
% completion supplies the prescribed six moments.
% opt.viscosity defaults to 1. Applied loads are divided by viscosity for
% the unit-viscosity solve; realisedForceTorque is returned in physical units.

if nargin == 0
    self_test();
    U = []; iters = []; lambda_norm = []; uerr = []; details = struct();
    return
end
if nargin < 3
    opt = struct();
end
opt = image_solver_options(opt);
validateattributes(q,{'numeric'},{'real','finite','2d','ncols',3});
P = size(q,1);
validateattributes(Fvec,{'numeric'}, ...
    {'real','finite','column','numel',6*P});
Fvec = Fvec/opt.viscosity;
hasPrepared = isfield(opt,'prepared') && ~isempty(opt.prepared);
if hasPrepared
    grid = validate_prepared_image_problem(q,opt,"mobility");
else
    grid = build_image_discretization(q,opt);
end
counts = arrayfun(@(body) size(body.collocation,1),grid.body);
counts = counts(:).';
offsets = [0,cumsum(3*counts)];
sourceCounts = arrayfun(@(body) 3*sum(arrayfun( ...
    @(b) size(b.points,1),body.basis)),grid.body);
sourceCounts = sourceCounts(:).';
sourceOffsets = [0,cumsum(sourceCounts)];
factors = repmat(struct('H',[],'C',[],'Y',[],'U',[], ...
    'rank',0,'singularValues',[],'sourceCount',0, ...
    'collocationCount',0),P,1);
H = cell(P,1);
C = cell(P,1);
lambda0 = cell(P,1);
coeff0 = cell(P,1);
for k = 1:P
    if opt.solverMode == "oneBody"
        if hasPrepared
            factors(k) = opt.prepared.mobilityFactors(k);
        else
            memoryMiB = 8*double(3*counts(k))*double(sourceCounts(k))/2^20;
            fprintf('mobility body %d: SVD %d x %d (%.1f MiB)\n', ...
                k,3*counts(k),sourceCounts(k),memoryMiB);
            factors(k) = one_body_image_mobility_preconditioner( ...
                grid.body(k),q(k,:),opt.svdCutoff);
        end
        H{k} = factors(k).H;
        C{k} = factors(k).C;
    else
        [H{k},C{k}] = image_mobility_maps(grid.body(k),q(k,:));
    end
    forceRows = 6*(k-1)+(1:6);
    lambda0{k} = C{k}*Fvec(forceRows);
    coeff0{k} = image_source_coefficients( ...
        grid.body(k).basis,lambda0{k});
end

rhs = zeros(offsets(end),1);
for k = 1:P
    target = grid.body(k);
    for first = 1:opt.batchSize:counts(k)
        nodes = first:min(first+opt.batchSize-1,counts(k));
        rows = offsets(k)+reshape(3*nodes-2+(0:2).',[],1);
        velocity = zeros(numel(nodes),3);
        for j = 1:P
            velocity = velocity+direct_image_velocity( ...
                grid.body(j).basis,coeff0{j},target.collocation(nodes,:));
        end
        rhs(rows) = -reshape((target.rowScale(nodes).*velocity).',[],1);
    end
end
z = cell(P,1);
if opt.solverMode == "oneBody"
    if norm(rhs) == 0
        mu = zeros(size(rhs));
        iters = 0;
        trueGmresResidual = 0;
    else
        fprintf('starting mobility GMRES\n');
        [mu,iters,~,~] = helsing_gmres( ...
            @(x) image_mobility_matvec(x,grid,factors,opt.batchSize), ...
            rhs,numel(rhs),opt.maxIterations,opt.gmresTolerance,0,0,true);
        fprintf('completed mobility GMRES\n');
        trueGmresResidual = norm(image_mobility_matvec( ...
            mu,grid,factors,opt.batchSize)-rhs)/norm(rhs);
    end
    for k = 1:P
        rows = offsets(k)+(1:3*counts(k));
        z{k} = factors(k).Y*(factors(k).U'*mu(rows));
    end
    retainedRanks = [factors.rank].';
    linearSystemResidual = trueGmresResidual;
else
    % Assemble the global recompleted system only for this dense reference.
    nRow = offsets(end);
    nCol = sourceOffsets(end);
    fprintf(['Dense mobility SVD matrix: %d rows x %d columns ', ...
        '(%.1f MiB for one double copy)\n'], ...
        nRow,nCol,8*double(nRow)*double(nCol)/2^20);
    matrix = zeros(nRow,nCol);
    for k = 1:P
        target = grid.body(k);
        for j = 1:P
            cols = sourceOffsets(j)+(1:sourceCounts(j));
            for first = 1:opt.batchSize:counts(k)
                nodes = first:min(first+opt.batchSize-1,counts(k));
                rows = offsets(k)+reshape(3*nodes-2+(0:2).',[],1);
                G = stokes_mfs_velocity_matrix( ...
                    grid.body(j).basis,target.collocation(nodes,:));
                block = G-(G*C{j})*H{j};
                if k == j
                    K = getKmat(target.collocation(nodes,:),q(k,:));
                    block = block+K*H{j};
                end
                matrix(rows,cols) = ...
                    repelem(target.rowScale(nodes),3).*block;
            end
        end
    end
    [Y,Ufactor,singularValues] = getPseudoFactors( ...
        matrix,opt.svdCutoff,0);
    assert(~isempty(singularValues),'The SVD cutoff removed every mode.');
    globalZ = Y*(Ufactor'*rhs);
    linearSystemResidual = norm(matrix*globalZ-rhs)/max(norm(rhs),eps);
    for k = 1:P
        cols = sourceOffsets(k)+(1:sourceCounts(k));
        z{k} = globalZ(cols);
    end
    retainedRanks = numel(singularValues);
    trueGmresResidual = NaN;
    iters = 0;
    clear matrix Y Ufactor singularValues globalZ
end

U = zeros(6*P,1);
lambda = cell(P,1);
realisedForceTorque = zeros(6*P,1);
lambda_norm = 0;
for k = 1:P
    forceRows = 6*(k-1)+(1:6);
    U(forceRows) = -H{k}*z{k};
    lambda{k} = z{k}-C{k}*(H{k}*z{k})+lambda0{k};
    realisedForceTorque(forceRows) = H{k}*lambda{k};
    lambda_norm = max(lambda_norm,max(abs(lambda{k})));
end
[uerr,validationPoints,validationVelocity,validationRelative, ...
    validationGap] = image_surface_validation(q,grid,lambda,U,opt);
details = struct('grid',grid,'lambda',{lambda},'z',{z}, ...
    'completion',{lambda0},'realisedForceTorque',opt.viscosity*realisedForceTorque, ...
    'retainedRanks',retainedRanks, ...
    'trueGmresResidual',trueGmresResidual, ...
    'linearSystemResidual',linearSystemResidual, ...
    'systemSize',[offsets(end),sourceOffsets(end)], ...
    'validationPoints',{validationPoints}, ...
    'validationVelocity',{validationVelocity}, ...
    'validationRelative',{validationRelative}, ...
    'validationGap',{validationGap},'options',opt);
end

function self_test()
%SELF_TEST checks moments, recovered motion, and an independent gap grid.
root = fileparts(fileparts(mfilename('fullpath')));
oldPath = path; cleanup = onCleanup(@() path(oldPath));
addpath(genpath(root));
q = [0,0,0;2.08,0,0];
Uref = [1;0.2;0;0;0;0.3;-0.7;0;0.1;0;0.2;0];
opt = struct('N',700,'collocationRatio',2,'N_image',4, ...
    'capLayers',[5,5],'capMinAzimuth',4, ...
    'validationCount',100,'validationGapCount',61, ...
    'gmresTolerance',1e-8,'maxIterations',150,'svdCutoff',1e-11,'viscosity',1);
[F,~,~,~,~] = solve_resistance_images(q,Uref,opt);
[U,~,~,err,details] = solve_mobility_images(q,F,opt);
assert(norm(details.realisedForceTorque-F,inf) < ...
    1e-10*max(1,norm(F,inf)));
assert(norm(U-Uref,inf) < 1e-3*max(1,norm(Uref,inf)));
assert(err < 0.1 && any(details.validationGap{1}));
for k = 1:2
    [H,C] = image_mobility_maps(details.grid.body(k),q(k,:));
    rows = 6*(k-1)+(1:6);
    assert(norm(U(rows)+H*details.z{k}) < 1e-11);
    assert(norm(H*details.completion{k}-F(rows)) < 1e-10);
    assert(norm(H*(details.z{k}-C*(H*details.z{k}))) < 1e-10);
end
basic = image_solver_options(struct('N',30,'N_image',0, ...
    'useImages',false,'collocationRatio',2));
grid = build_image_discretization(q,basic);
assert(all(arrayfun(@(body) isempty(body.capPoints),grid.body)));
[H,~] = image_mobility_maps(grid.body(1),q(1,:));
assert(norm(H-getKmat(grid.body(1).basis.points,q(1,:)).') == 0);
q3 = [0,0,0;2.08,0,0;4.16,0,0];
small = opt;
small.N = 30;
small.N_image = 2;
small.useRowScaling = false;
F3 = sin((1:18).');
[U3,~,~,~,d3] = solve_mobility_images(q3,F3,small);
assert(all(isfinite(U3)) && ...
    norm(d3.realisedForceTorque-F3,inf) < 1e-10);
assert(all(arrayfun(@(body) all(body.rowScale == 1),d3.grid.body)));
assert(size(d3.grid.body(2).collocation,1) > ...
    size(d3.grid.body(1).collocation,1));
assert(numel(d3.lambda{2}) > numel(d3.lambda{1}));
dense = small;
dense.solverMode = "denseSVD";
[Ud,itd,~,~,dd] = solve_mobility_images(q,F,dense);
assert(itd == 0 && isnan(dd.trueGmresResidual));
assert(isscalar(dd.retainedRanks) && dd.retainedRanks > 0);
assert(norm(dd.realisedForceTorque-F,inf) < 1e-10);
weightedError = cell(2,1);
weightedRhs = cell(2,1);
for k = 1:2
    body = dd.grid.body(k);
    v = zeros(3*size(body.collocation,1),1);
    v0 = v;
    for j = 1:2
        G = stokes_mfs_velocity_matrix( ...
            dd.grid.body(j).basis,body.collocation);
        v = v+G*dd.lambda{j};
        v0 = v0+G*dd.completion{j};
    end
    rows = 6*(k-1)+(1:6);
    w = repelem(body.rowScale,3);
    weightedError{k} = w.*(v-getKmat(body.collocation,q(k,:))*Ud(rows));
    weightedRhs{k} = -w.*v0;
end
physicalResidual = norm(vertcat(weightedError{:}))/ ...
    norm(vertcat(weightedRhs{:}));
assert(abs(physicalResidual-dd.linearSystemResidual) < 1e-11);
fprintf('solve_mobility_images self-test passed: err %.3e.\n',err);
end
