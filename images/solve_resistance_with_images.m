function [Fvec,iters,lambda_norm,err_res,details] = ...
    solve_resistance_with_images(q,U,opt)
%SOLVE_RESISTANCE_WITH_IMAGES Solve rigid-sphere Stokes resistance by image MFS.
%   [Fvec,iters,lambda_norm,err_res,details] = ...
%       solve_resistance_with_images(q,U,opt)
%   solve_resistance_with_images() shows a three-sphere resistance example.
% Reference: A. Broms, A. H. Barnett, and A.-K. Tornberg,
% "Accurate close interactions of Stokes spheres using lubrication-adapted
% image systems," Journal of Computational Physics 523 (2025), 113636.
%
% Inputs: q is P-by-3 of unit-sphere centres; U is 6P-by-1 with each body
% [vx;vy;vz;omega_x;omega_y;omega_z]. opt is a small flat control struct;
% omitted fields receive the shared image-solver defaults. Image families can be
% any nonempty subset of SRTD; coarse proxy nodes always carry S.
% Main controls: N (coarse proxy count), proxyRadius, nearThreshold,
% N_image (nodes per close-contact image line; [] uses paper Eq. 39),
% imageFamilies, capMode, capAngles, capUpsampling or exact capPointCounts,
% capLayers, capExponent, rowWeightMode, useRowScaling (default true),
% svdCutoff, gmresTolerance and
% solverMode ("oneBody" by default, or "denseSVD").
% Outputs: Fvec is 6P-by-1 [Fx;Fy;Fz;Tx;Ty;Tz] per body, using the same
% source-strength sign as solve_resistance/Kin'*lambda. iters is the GMRES
% count (zero for denseSVD), lambda_norm the maximum absolute source
% strength, and err_res the
% maximum independent pointwise relative velocity error. At validation
% nodes with zero prescribed velocity, the largest prescribed speed is
% used as the denominator. details contains
% the grids, normalized source coefficients, retained ranks, true
% linear-system residual,
% and independent validation locations, computed velocities, and relative
% residuals. For the
% dense solve, retainedRanks is the global truncated rank and
% trueGmresResidual is NaN; linearSystemResidual is the true weighted
% least-squares residual.
%
% For each sphere, B_k = W_k*G_kk is SVD factored with relative truncation.
% GMRES solves an identity-diagonal right-preconditioned system in weighted
% collocation velocity unknowns. Off-diagonal interactions are direct sums
% in batches. This default path does not assemble the global matrix. The
% optional denseSVD path assembles the same complete weighted collocation
% system and applies its truncated pseudoinverse through two SVD factors.
% During one-body setup, progress reports give each body's SVD dimensions
% and the storage needed for one double-precision copy of that matrix.
% The velocity solve uses unit-viscosity coefficients. For physical viscosity
% opt.viscosity, S/R/T strengths and the returned force/torque scale by it;
% D strengths do not. details.lambda and lambda_norm remain normalized
% unit-viscosity coefficients, suitable for direct_image_velocity.
% The sphere radius is one.

if nargin == 0
    self_test();
    Fvec = []; iters = []; lambda_norm = []; err_res = []; details = struct();
    return
end
if nargin < 3
    opt = struct();
end
opt = image_solver_options(opt);
validateattributes(q,{'numeric'},{'real','finite','2d','ncols',3});
P = size(q,1);
validateattributes(U,{'numeric'},{'real','finite','column','numel',6*P});

grid = build_image_discretization(q,opt);
rhsBlocks = cell(P,1);
for k = 1:P
    body = grid.body(k);
    physicalData = getKmat(body.collocation,q(k,:))*U(6*(k-1)+(1:6));
    rhsBlocks{k} = repelem(body.rowScale,3).*physicalData;
end
rhs = vertcat(rhsBlocks{:});
lambda = cell(P,1);
collocationCounts = arrayfun( ...
    @(body) size(body.collocation,1),grid.body);
collocationCounts = collocationCounts(:).';
collocationOffsets = [0,cumsum(3*collocationCounts)];
N_total = arrayfun(@(body) sum(arrayfun( ...
    @(block) size(block.points,1),body.basis)),grid.body);
N_total = N_total(:).';
N_offsets = [0,cumsum(3*N_total)];
if opt.solverMode == "oneBody"
    factors = repmat(struct('Y',[],'U',[],'rank',0, ...
        'singularValues',[],'N_total',0,'collocationCount',0),P,1);
    for k = 1:P
        if k > 1 && ~opt.useImages && ...
                isequal(size(grid.body(k).collocation), ...
                    size(grid.body(1).collocation)) && ...
                isscalar(grid.body(k).basis) && ...
                isscalar(grid.body(1).basis) && ...
                grid.body(k).basis(1).family == "S" && ...
                grid.body(1).basis(1).family == "S" && ...
                grid.body(k).basis(1).scale == ...
                    grid.body(1).basis(1).scale && ...
                isequal(grid.body(k).rowScale,grid.body(1).rowScale) && ...
                max(abs(grid.body(k).collocation-q(k,:) - ...
                    (grid.body(1).collocation-q(1,:))),[],'all') < 1e-14 && ...
                max(abs(grid.body(k).basis(1).points-q(k,:) - ...
                    (grid.body(1).basis(1).points-q(1,:))),[],'all') < 1e-14
            % Translated identical proxy spheres have the same self block.
            factors(k) = factors(1);
            fprintf('reuse preconditioner for body %d from body 1\n',k);
            continue
        end
        nRow = 3*collocationCounts(k);
        nCol = 3*N_total(k);
        memMiB = 8*double(nRow)*double(nCol)/2^20;
        fprintf(['determine preconditioner for body %d: SVD matrix ', ...
            '%d rows x %d cols (%.1f MiB, double storage)\n'], ...
            k,nRow,nCol,memMiB);
        factors(k) = one_body_image_preconditioner( ...
            grid.body(k),opt.svdCutoff);
    end
    if norm(rhs) == 0
        mu = zeros(size(rhs));
        iters = 0;
        trueGmresResidual = 0;
    else
        fprintf('starting GMRES\n');
        [mu,iters,~,~] = helsing_gmres( ...
            @(x) image_resistance_matvec(x,grid,factors,opt.batchSize), ...
            rhs,numel(rhs),opt.maxIterations,opt.gmresTolerance,0,0,true);
        fprintf('completed GMRES\n');
        trueGmresResidual = norm(image_resistance_matvec( ...
            mu,grid,factors,opt.batchSize)-rhs)/norm(rhs);
    end
    for k = 1:P
        rows = collocationOffsets(k)+(1:3*collocationCounts(k));
        lambda{k} = factors(k).Y*(factors(k).U'*mu(rows));
    end
    retainedRanks = [factors.rank].';
    linearSystemResidual = trueGmresResidual;
else
    % Each target/source pair fills one weighted block of the global
    % rectangular matrix. Target batching bounds temporary kernel storage.
    nRow = collocationOffsets(end);
    nCol = N_offsets(end);
    fprintf(['Dense resistance SVD matrix: %d rows x %d columns ', ...
        '(%.1f MiB for one double copy)\n'], ...
        nRow,nCol,8*double(nRow)*double(nCol)/2^20);
    matrix = zeros(collocationOffsets(end),N_offsets(end));
    for targetBody = 1:P
        target = grid.body(targetBody);
        for sourceBody = 1:P
            columns = N_offsets(sourceBody)+(1:3*N_total(sourceBody));
            for first = 1:opt.batchSize:collocationCounts(targetBody)
                nodes = first:min(first+opt.batchSize-1, ...
                    collocationCounts(targetBody));
                rows = collocationOffsets(targetBody)+ ...
                    reshape(3*nodes-2+(0:2).',[],1);
                kernel = stokes_mfs_velocity_matrix( ...
                    grid.body(sourceBody).basis,target.collocation(nodes,:));
                matrix(rows,columns) = ...
                    repelem(target.rowScale(nodes),3).*kernel;
            end
        end
    end
    [Y,Ufactor,singularValues] = getPseudoFactors( ...
        matrix,opt.svdCutoff,0);
    assert(~isempty(singularValues),'The SVD cutoff removed every mode.');
    globalStrength = Y*(Ufactor'*rhs);
    linearSystemResidual = norm(matrix*globalStrength-rhs)/max(norm(rhs),eps);
    for k = 1:P
        columns = N_offsets(k)+(1:3*N_total(k));
        lambda{k} = globalStrength(columns);
    end
    retainedRanks = numel(singularValues);
    trueGmresResidual = NaN;
    iters = 0;
    clear matrix Y Ufactor singularValues globalStrength
end

Fvec = zeros(6*P,1);
lambda_norm = 0;
for k = 1:P
    if ~isempty(lambda{k})
        lambda_norm = max(lambda_norm,max(abs(lambda{k})));
    end
    offset = 0;
    for blockIndex = 1:numel(grid.body(k).basis)
        block = grid.body(k).basis(blockIndex);
        count = size(block.points,1);
        columns = offset+(1:3*count);
        strength = block.scale*reshape(lambda{k}(columns),3,[]).';
        if block.family == "S"
            Fvec(6*(k-1)+(1:3)) = Fvec(6*(k-1)+(1:3))+ ...
                sum(strength,1).';
            Fvec(6*(k-1)+(4:6)) = Fvec(6*(k-1)+(4:6))+ ...
                sum(cross(block.points-q(k,:),strength,2),1).';
        elseif block.family == "R"
            Fvec(6*(k-1)+(4:6)) = Fvec(6*(k-1)+(4:6))+ ...
                sum(strength,1).';
        end
        offset = offset+3*count;
    end
end
Fvec = opt.viscosity*Fvec;

% Use the same independent full-sphere and gap-cap check as mobility.
[err_res,validationPoints,validationVelocity,validationRelative, ...
    validationGap] = image_surface_validation(q,grid,lambda,U,opt);
details = struct('grid',grid,'lambda',{lambda}, ...
    'retainedRanks',retainedRanks, ...
    'trueGmresResidual',trueGmresResidual, ...
    'linearSystemResidual',linearSystemResidual, ...
    'systemSize',[collocationOffsets(end),N_offsets(end)], ...
    'validationPoints',{validationPoints}, ...
    'validationVelocity',{validationVelocity}, ...
    'validationRelative',{validationRelative}, ...
    'validationGap',{validationGap},'options',opt);
end

function self_test()
%SELF_TEST displays paper-scale three-sphere velocities and resistance data.
fprintf(['Three-sphere resistance self-test: prescribe inward rigid ', ...
    'translations for a close triangle, then solve for forces, torques, ', ...
    'and surface velocities with gap-adaptive image lines.\n']);
fprintf(['Reference: Broms, Barnett and Tornberg, "Accurate close ', ...
    'interactions of Stokes spheres using lubrication-adapted image ', ...
    'systems" (Journal of Computational Physics 523, 2025).\n']);
root = fileparts(fileparts(mfilename('fullpath')));
oldPath = path; cleanup = onCleanup(@() path(oldPath));
addpath(genpath(root));
gap = 1e-3;
side = 2+gap;
q = [0,0,0;side,0,0;side/2,sqrt(3)*side/2,0];
U = zeros(18,1);
for k = 1:3
    direction = mean(q,1)-q(k,:);
    U(6*(k-1)+(1:3)) = (direction/norm(direction)).';
end
N = 686;
N_image = []; % Adaptive from each gap by paper Eq. (39).
opt = struct('N',N,'proxyRadius',0.63,'N_image',N_image, ...
    'capMode',"exponential",'rowWeightMode',"local", ...
    'validationCount',1000,'validationGapCount',401);
[F,iters,lambdaNorm,residual,details] = ...
    solve_resistance_with_images(q,U,opt);
fprintf('Three-sphere triangle, gap %.1e, N=%d, N_image=%d per pair\n', ...
    gap,details.grid.N,details.grid.pairs(1,4));
fprintf('GMRES iterations: %d; true weighted residual: %.3e\n', ...
    iters,details.trueGmresResidual);
fprintf('Maximum independent relative surface residual: %.3e\n',residual);
fprintf('Maximum source strength: %.3e\n',lambdaNorm);
for k = 1:size(q,1)
    rows = 6*(k-1)+(1:6);
    gap = details.validationGap{k};
    fprintf(['Body %d: F = [% .4g % .4g % .4g], ', ...
        'T = [% .4g % .4g % .4g], surface %.3e, gap %.3e\n'], ...
        k,F(rows),max(details.validationRelative{k}), ...
        max(details.validationRelative{k}(gap)));
end

fig = figure('Name','Three-sphere image resistance: surface velocity');
layout = tiledlayout(fig,1,2,'TileSpacing','compact');
speedAx = nexttile(layout);
resAx = nexttile(layout);
hold(speedAx,'on');
hold(resAx,'on');
for k = 1:size(q,1)
    points = details.validationPoints{k};
    velocity = details.validationVelocity{k};
    speed = vecnorm(velocity,2,2);
    scatter3(speedAx,points(:,1),points(:,2),points(:,3), ...
        11,speed,'filled');
    % Arrows use only the independent full-sphere nodes for legibility.
    nFull = nnz(~details.validationGap{k});
    sample = 1:max(1,floor(nFull/28)):nFull;
    quiver3(speedAx,points(sample,1),points(sample,2),points(sample,3), ...
        velocity(sample,1),velocity(sample,2),velocity(sample,3), ...
        0.35,'Color',[0.1,0.1,0.1]);
    scatter3(resAx,points(:,1),points(:,2),points(:,3), ...
        11,log10(max(details.validationRelative{k},eps)),'filled');
end
for ax = [speedAx,resAx]
    axis(ax,'equal');
    grid(ax,'on');
    view(ax,3);
    xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
    colorbar(ax);
end
title(speedAx,'Computed surface speed and velocity direction');
title(resAx,'log_{10} relative boundary residual');
end
