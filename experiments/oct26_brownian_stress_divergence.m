function results = oct26_brownian_stress_divergence()
%OKT26_BROWNIAN_STRESS_DIVERGENCE Two-sphere image-free MFS stresslet sweep.
% Run deliberately: the complete sweep, refinements and reference solves
% are substantial. The plotted quantity is D, not k_B*T*D. With the
% force-on-fluid convention the deviatoric suspension stress is
% -(k_B*T/V)*sum_p D_p.

%% User controls
% Geometry and discretization.
radius = 1;
viscosity = 1;                 % Basic solvers use unit viscosity.
requestedSources = 686;
spatialSources = [302,686,1000];
proxyRadius = 0.68;
collocationRatio = 1.2;
gapRatios = logspace(-2,0,9);

% Solver tolerances and finite-difference settings.
gmresTolerance = 1e-8;
maxIterations = 400;
surfaceTolerance = 1e-3;
stressletTolerance = 1e-2;
divergenceTolerance = 0.05;
angleStep = 0.01;
stepLevels = 4;

% Optional validation work and representative gaps.
runRotationalStudy = true;
runFullMatrixReference = true;
runSpatialCheckAtEveryGap = true;
studyIndices = [1,5,9];       % Wide, intermediate and narrow.

assert(radius==1 && viscosity==1, ...
    'The selected basic sphere solvers use radius and viscosity one.');
settings = struct('radius',radius,'viscosity',viscosity, ...
    'requestedSources',requestedSources,'spatialSources',spatialSources, ...
    'proxyRadius',proxyRadius,'collocationRatio',collocationRatio, ...
    'gapRatios',gapRatios,'gmresTolerance',gmresTolerance, ...
    'maxIterations',maxIterations, ...
    'surfaceTolerance',surfaceTolerance, ...
    'stressletTolerance',stressletTolerance, ...
    'divergenceTolerance',divergenceTolerance, ...
    'angleStep',angleStep,'stepLevels',stepLevels, ...
    'runRotationalStudy',runRotationalStudy, ...
    'runFullMatrixReference',runFullMatrixReference, ...
    'runSpatialCheckAtEveryGap',runSpatialCheckAtEveryGap, ...
    'studyIndices',studyIndices,'matlabVersion',version);
root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(root,'setpath.m'));
opt = init_MFS(requestedSources);
opt.a_glob = collocationRatio;
opt.ellipsoid = true;          % Body-frame path also handles rotated spheres.
opt.lr = 0;
opt.fmm = 0;
opt.plot = false;
opt.debug = false;
opt.profile = false;
opt.gmres_verbose = 0;
opt.compute_residual = true;
opt.gmres_tol = gmresTolerance;
opt.maxit = maxIterations;

%% One-body work, once per spatial resolution
% Cache body-frame grids and one-body preconditioners for each unique grid.
levels = unique(spatialSources,'stable');
setups = cell(numel(levels),1);
for level = 1:numel(levels)
    gridOpt = opt;
    gridOpt.des_n = levels(level);
    [rin0,rout0] = getDesignGrid(proxyRadius,gridOpt);
    setup = struct();
    setup.rin0 = rin0;
    setup.rout0 = rout0;
    setup.cache = prepare_sphere_mfs_preconditioners(rin0,rout0);
    setup.requestedSources = levels(level);
    setups{level} = setup;
    fprintf('N requested %d, actual %d, M %d: setup %.2f s\n', ...
        levels(level),size(rin0,1),size(rout0,1), ...
        setup.cache.setupSeconds);
end
baseIndex = find(levels==requestedSources,1);
assert(~isempty(baseIndex),'requestedSources must occur in spatialSources.');
base = setups{baseIndex};
setupSummary = repmat(struct('requestedSources',0,'actualSources',0, ...
    'collocationCount',0,'seconds',0,'resistanceRank',0, ...
    'mobilityRank',0),numel(levels),1);
for level = 1:numel(levels)
    setup = setups{level};
    setupSummary(level) = struct( ...
        'requestedSources',setup.requestedSources, ...
        'actualSources',size(setup.rin0,1), ...
        'collocationCount',size(setup.rout0,1), ...
        'seconds',setup.cache.setupSeconds, ...
        'resistanceRank',size(setup.cache.resistance.Y,2), ...
        'mobilityRank',size(setup.cache.mobility.Y,2));
end
nGaps = numel(gapRatios);
results = cell(nGaps,1);
stepStudies = cell(numel(studyIndices),1);
selectedSteps = NaN(numel(studyIndices),1);

%% Establish finite-difference steps before the production sweep
% At representative gaps, halve the translation step and select the first
% finer level whose change is below tolerance; fall back to the finest level.
for t = 1:numel(studyIndices)
    k = studyIndices(t);
    gap = gapRatios(k)*radius;
    q = pair_centres(radius,gap);
    Q = {eye(3),eye(3)};
    h0 = min(0.01*radius,0.04*gap);
    study = repmat(struct('h',0,'D',[],'Dbody',[], ...
        'residual',0,'seconds',0),stepLevels,1);
    for level = 1:stepLevels
        h = h0/2^(level-1);
        t0 = tic;
        div = centered_divergence(q,Q,base,opt,h,angleStep,false);
        study(level).h = h;
        study(level).D = sum(div.translation,3);
        study(level).Dbody = div.translation;
        study(level).residual = div.maxSurfaceResidual;
        study(level).seconds = toc(t0);
    end
    changes = zeros(stepLevels-1,1);
    for level = 1:stepLevels-1
        difference = study(level).Dbody-study(level+1).Dbody;
        scale = max(norm(study(level+1).Dbody(:)),1e-12);
        changes(level) = norm(difference(:))/scale;
    end
    observedRatios = changes(1:end-1)./max(changes(2:end),eps);
    good = find(changes<=divergenceTolerance,1,'first');
    if isempty(good)
        selectedSteps(t) = study(end).h;
    else
        selectedSteps(t) = study(good+1).h;
    end
    stepStudies{t} = struct('gap',gap,'levels',study, ...
        'successiveChanges',changes,'observedRatios',observedRatios, ...
        'selectedStep',selectedSteps(t));
end

%% Main sweep: raw maps and complete-product finite differences
% At each gap, measure the maps, estimate step and spatial errors, apply the
% acceptance checks, then save a checkpoint so completed gaps are retained.
for k = 1:nGaps
    gap = gapRatios(k)*radius;
    q = pair_centres(radius,gap);
    Q = {eye(3),eye(3)};
    [~,nearest] = min(abs(log(gapRatios(k))-log(gapRatios(studyIndices))));
    h = min(selectedSteps(nearest),0.04*gap);
    entry = struct('gap',gap,'gapOverRadius',gap/radius, ...
        'translationStep',h,'angleStep',angleStep, ...
        'requestedSources',requestedSources,'actualSources',size(base.rin0,1), ...
        'collocationCount',size(base.rout0,1),'base',[], ...
        'divergence',[],'stepError',NaN,'spatialError',NaN, ...
        'rotation',[],'fullMatrixReference',[], ...
        'valid',false,'failure','');
    try
        entry.base = measure_maps(q,Q,base,opt);
        entry.divergence = centered_divergence( ...
            q,Q,base,opt,h,angleStep,false);
        Dbody = entry.divergence.translation;
        Dsum = sum(Dbody,3);
        offDiagonal = Dsum-diag(diag(Dsum));
        entry.symmetryDiagnostics = struct( ...
            'yyMinusZz',Dsum(2,2)-Dsum(3,3), ...
            'offDiagonalFrobenius',norm(offDiagonal,'fro'), ...
            'trace',trace(Dsum), ...
            'normalizer',max(norm(Dsum,'fro'),1e-12));
        halfStep = centered_divergence( ...
            q,Q,base,opt,h/2,angleStep,false);
        D = entry.divergence.translation;
        Dhalf = halfStep.translation;
        stepDifference = D-Dhalf;
        entry.stepDifference = norm(stepDifference(:));
        entry.stepHalf = halfStep;
        tensorScale = max([norm(D(:)),norm(Dhalf(:)),1e-12]);
        if runSpatialCheckAtEveryGap || ismember(k,studyIndices)
            fineIndex = find(levels==max(spatialSources),1);
            fine = centered_divergence( ...
                q,Q,setups{fineIndex},opt,h,angleStep,false);
            Dfine = fine.translation;
            spatialDifference = D-Dfine;
            entry.spatialDifference = norm(spatialDifference(:));
            tensorScale = max(tensorScale,norm(Dfine(:)));
            entry.spatial = fine;
            if ismember(k,studyIndices)
                coarseIndex = find(levels==min(spatialSources),1);
                coarse = centered_divergence( ...
                    q,Q,setups{coarseIndex},opt,h,angleStep,false);
                Dcoarse = coarse.translation;
                coarseDifference = Dcoarse-D;
                entry.coarseSpatialError = norm(coarseDifference(:))/ ...
                    max(norm(D(:)),1e-12);
                entry.spatialCoarse = coarse;
            end
        end
        % Use the same per-gap Frobenius tensor scale for both refinements.
        entry.refinementScale = tensorScale;
        entry.stepError = entry.stepDifference/tensorScale;
        if isfield(entry,'spatialDifference')
            entry.spatialError = entry.spatialDifference/tensorScale;
        end
        failedChecks = {};
        if ~all(isfinite([entry.base.maxSurfaceResidual, ...
                entry.divergence.maxSurfaceResidual, ...
                halfStep.maxSurfaceResidual, ...
                entry.base.stressletConsistency,entry.stepError]))
            failedChecks{end+1} = 'nonfinite diagnostic'; %#ok<AGROW>
        end
        if entry.base.maxSurfaceResidual>surfaceTolerance || ...
                entry.divergence.maxSurfaceResidual>surfaceTolerance || ...
                halfStep.maxSurfaceResidual>surfaceTolerance || ...
                (isfield(entry,'spatial') && ...
                entry.spatial.maxSurfaceResidual>surfaceTolerance)
            failedChecks{end+1} = 'surface residual'; %#ok<AGROW>
        end
        if entry.base.stressletConsistency>stressletTolerance
            failedChecks{end+1} = 'stresslet consistency'; %#ok<AGROW>
        end
        if any(entry.base.forceMomentError>1e-8) || ...
                any(~isfinite(entry.base.forceMomentError))
            failedChecks{end+1} = 'force/torque moment'; %#ok<AGROW>
        end
        if entry.stepError>divergenceTolerance
            failedChecks{end+1} = 'step refinement'; 
        end
        if ~isfinite(entry.spatialError) || ...
                entry.spatialError>divergenceTolerance
            failedChecks{end+1} = 'spatial refinement unavailable or failed';
        end
        entry.valid = isempty(failedChecks);
        entry.failure = strjoin(failedChecks,', ');
    catch ME
        entry.failure = sprintf('%s: %s',ME.identifier,ME.message);
        warning('Gap %.4g failed: %s',gap,entry.failure);
    end
    results{k} = entry;
    fprintf('gap/a %.4g: valid=%d, step %.3e, spatial %.3e\n', ...
        gap/radius,entry.valid,entry.stepError,entry.spatialError);
    save(fullfile(root,'experiments', ...
        'okt26_brownian_stress_divergence_results.mat'), ...
        'results','stepStudies','setupSummary','settings','gapRatios', ...
        'radius','viscosity','proxyRadius','collocationRatio','opt','-v7.3');
end

%% Selected full-product and orientation checks
% More expensive cross-checks run only at representative gaps that passed
% the production checks above.
for t = 1:numel(studyIndices)
    k = studyIndices(t);
    if ~results{k}.valid
        continue
    end
    q = pair_centres(radius,results{k}.gap);
    Q = {eye(3),eye(3)};
    h = results{k}.translationStep;
    if runFullMatrixReference
        reference = full_matrix_difference(q,Q,base,opt,h);
        direct = results{k}.divergence.translation;
        referenceDifference = direct-reference.Dbody;
        reference.relativeDifference = norm(referenceDifference(:))/ ...
            max(norm(direct(:)),1e-12);
        results{k}.fullMatrixReference = reference;
        results{k}.valid = results{k}.valid && ...
            reference.relativeDifference<=stressletTolerance && ...
            reference.maxSurfaceResidual<=surfaceTolerance;
        if ~results{k}.valid
            results{k}.failure = 'explicit full-product reference';
        end
    end
    if runRotationalStudy
        rotation = cell(numel(levels),1);
        for level = 1:numel(levels)
            rotation{level} = centered_divergence( ...
                q,Q,setups{level},opt,h,angleStep,true);
        end
        results{k}.rotation = rotation;
        DrotBase = sum(rotation{baseIndex}.rotation,3);
        DrotFine = sum(rotation{find(levels==max(spatialSources),1)}.rotation,3);
        results{k}.rotationSpatialError = norm(DrotBase-DrotFine,'fro')/ ...
            max(norm(DrotFine,'fro'),1e-12);
        halfAngle = centered_divergence( ...
            q,Q,base,opt,h,angleStep/2,true);
        results{k}.rotationHalfAngle = halfAngle;
        results{k}.rotationStepError = norm(DrotBase- ...
            sum(halfAngle.rotation,3),'fro')/ ...
            max(norm(sum(halfAngle.rotation,3),'fro'),1e-12);
        results{k}.fullD = results{k}.divergence.translation+ ...
            rotation{baseIndex}.rotation;
    end
end

%% Signed plot of the measured, summed translational divergence
% Failed gaps remain NaN and are omitted from the plotted curves.
figureHandle = figure('Name','Two-sphere stresslet divergence');
ax = axes(figureHandle);
hold(ax,'on'); grid(ax,'on');
diagonal = NaN(nGaps,3);
for k = 1:nGaps
    if results{k}.valid
        D = sum(results{k}.divergence.translation,3);
        diagonal(k,:) = diag(D).';
    end
end
plot(ax,gapRatios,diagonal(:,1),'o-','DisplayName','D_{xx}');
plot(ax,gapRatios,diagonal(:,2),'s-','DisplayName','D_{yy}');
plot(ax,gapRatios,diagonal(:,3),'^-','DisplayName','D_{zz}');
set(ax,'XScale','log');
xlabel(ax,'Surface gap \delta/a');
ylabel(ax,'Summed D (force-on-fluid moment; no k_B T)');
legend(ax,'Location','best');
title(ax,'Image-free MFS: translational stresslet divergence');
exportgraphics(figureHandle,fullfile(root,'visuals', ...
    'okt26_brownian_stress_divergence.png'));
save(fullfile(root,'experiments', ...
    'okt26_brownian_stress_divergence_results.mat'), ...
    'results','stepStudies','setupSummary','settings','gapRatios', ...
    'radius','viscosity','proxyRadius','collocationRatio','opt','-v7.3');
end

function q = pair_centres(radius,gap)
%PAIR_CENTRES places two equal spheres symmetrically on the x axis.
distance = radius+gap/2;
q = [-distance,0,0;distance,0,0];
end

function [rin,rout] = placed_grids(q,Q,setup)
%PLACED_GRIDS moves fixed body-frame nodes without remeshing.
P = size(q,1);
N = size(setup.rin0,1);
M = size(setup.rout0,1);
rin = zeros(N*P,3);
rout = zeros(M*P,3);
for p = 1:P
    rin((p-1)*N+(1:N),:) = (Q{p}*setup.rin0.').'+q(p,:);
    rout((p-1)*M+(1:M),:) = (Q{p}*setup.rout0.').'+q(p,:);
end
end

function maps = measure_maps(q,Q,setup,opt)
%MEASURE_MAPS records the resistance stresslet map and raw body mobility.
[rin,rout] = placed_grids(q,Q,setup);
P = size(q,1);
d = 6*P;
Rsu = zeros(5*P,d);
Rtensor = zeros(3,3,P,d);
Mbody = zeros(d,d);
mobilityStress = zeros(5*P,d);
mobilityStressTensor = zeros(3,3,P,d);
residualResistance = zeros(d,1);
residualMobility = zeros(d,1);
iterationsResistance = zeros(d,1);
iterationsMobility = zeros(d,1);
timeResistance = zeros(d,1);
timeMobility = zeros(d,1);
forceMomentError = zeros(d,1);
optRes = opt; optRes.precond = setup.cache.resistance;
optMob = opt; optMob.precond = setup.cache.mobility;
for j = 1:d
    e = zeros(d,1); e(j) = 1;
    timer = tic;
    [F,it,~,err,lambda] = solve_resistance( ...
        q,rin,rout,e,optRes,Q,[1,1,1]);
    timeResistance(j) = toc(timer);
    [S,s5] = apply_KS_transpose(rin,q,lambda);
    Rsu(:,j) = s5(:);
    Rtensor(:,:,:,j) = S;
    residualResistance(j) = err;
    iterationsResistance(j) = it;
    % For the deterministic resistance solve, moments equal reported loads.
    N = size(setup.rin0,1);
    for p = 1:P
        nodes = (p-1)*N+(1:N);
        rows = 3*(nodes(1)-1)+1:3*nodes(end);
        moment = getKmat(rin(nodes,:),q(p,:)).'*lambda(rows);
        assert(norm(moment-F(6*(p-1)+(1:6)))< ...
            1e-9*max(1,norm(F)));
    end
end
for j = 1:d
    e = zeros(d,1); e(j) = 1;
    timer = tic;
    [U,it,~,err,lambda] = solve_mobility( ...
        q,rin,rout,e,[],optMob,Q,[1,1,1]);
    timeMobility(j) = toc(timer);
    [Sdirect,s5] = apply_KS_transpose(rin,q,lambda);
    Mbody(:,j) = U;
    mobilityStress(:,j) = s5(:);
    mobilityStressTensor(:,:,:,j) = Sdirect;
    residualMobility(j) = err;
    iterationsMobility(j) = it;
    N = size(setup.rin0,1);
    for p = 1:P
        nodes = (p-1)*N+(1:N);
        rows = 3*(nodes(1)-1)+1:3*nodes(end);
        recovered = getKmat(rin(nodes,:),q(p,:)).'*lambda(rows);
        localForce = e(6*(p-1)+(1:6));
        forceMomentError(j) = max(forceMomentError(j), ...
            norm(recovered-localForce)/max(norm(localForce),1));
    end
end
A = Rsu*Mbody;
Atensor = reshape(reshape(Rtensor,9*P,d)*Mbody,3,3,P,d);
maps = struct('R_SU',Rsu,'R_SU_tensor',Rtensor, ...
    'Mbody',Mbody,'A',A,'A_tensor',Atensor, ...
    'mobilityStress',mobilityStress, ...
    'mobilityStressTensor',mobilityStressTensor, ...
    'stressletConsistency',norm(mobilityStress-A,'fro')/ ...
        max(norm(mobilityStress,'fro'),1e-12), ...
    'rawMobilitySymmetryDefect',norm(Mbody-Mbody.','fro')/ ...
        max(norm(Mbody,'fro'),1e-12), ...
    'forceMomentError',forceMomentError, ...
    'residualResistance',residualResistance, ...
    'residualMobility',residualMobility, ...
    'maxSurfaceResidual',max([residualResistance;residualMobility]), ...
    'iterationsResistance',iterationsResistance, ...
    'iterationsMobility',iterationsMobility, ...
    'timeResistance',timeResistance,'timeMobility',timeMobility, ...
    'rightHandSidesResistance',d,'rightHandSidesMobility',d);
end

function div = centered_divergence(q,Q,setup,opt,h,hAngle,includeRotation)
%CENTERED_DIVERGENCE differences stresslets from unit-force mobility solves.
% Translation and optional rotation perturbations are evaluated at -step/+step.
P = size(q,1);
d = 6*P;
Dtrans = zeros(3,3,P);
Drot = zeros(3,3,P);
residuals = zeros(2*d,1);
iterations = zeros(2*d,1);
solveTimes = zeros(2*d,1);
count = 0;
for j = 1:d
    component = mod(j-1,6)+1;
    if component>3 && ~includeRotation
        continue
    end
    body = floor((j-1)/6)+1;
    step = h;
    if component>3
        step = hAngle;
    end
    stressMinus = [];
    stressPlus = [];
    for perturbationSign = [-1,1]
        shiftedQ = q;
        shiftedR = Q;
        if component<=3
            shiftedQ(body,component) = shiftedQ(body,component)+ ...
                perturbationSign*step;
        else
            axis = zeros(3,1); axis(component-3) = 1;
            crossMatrix = [0,-axis(3),axis(2); ...
                axis(3),0,-axis(1);-axis(2),axis(1),0];
            shiftedR{body} = expm(perturbationSign*step*crossMatrix)*Q{body};
        end
        assert(norm(shiftedQ(2,:)-shiftedQ(1,:))>2, ...
            'Perturbation would overlap the spheres.');
        [rin,rout] = placed_grids(shiftedQ,shiftedR,setup);
        F = zeros(d,1); F(j) = 1;
        optMob = opt; optMob.precond = setup.cache.mobility;
        timer = tic;
        [~,it,~,err,lambda] = solve_mobility( ...
            shiftedQ,rin,rout,F,[],optMob,shiftedR,[1,1,1]);
        count = count+1;
        solveTimes(count) = toc(timer);
        iterations(count) = it;
        residuals(count) = err;
        stress = apply_KS_transpose(rin,shiftedQ,lambda);
        if perturbationSign < 0
            stressMinus = stress;
        else
            stressPlus = stress;
        end
    end
    derivative = (stressPlus-stressMinus)/(2*step);
    if component<=3
        Dtrans = Dtrans+derivative;
    else
        Drot = Drot+derivative;
    end
end
div = struct('translation',Dtrans,'rotation',Drot, ...
    'sumTranslation',sum(Dtrans,3),'sumRotation',sum(Drot,3), ...
    'residuals',residuals(1:count), ...
    'maxSurfaceResidual',max(residuals(1:count)), ...
    'iterations',iterations(1:count), ...
    'solveTimes',solveTimes(1:count), ...
    'mobilityRightHandSides',count,'stepTranslation',h, ...
    'stepRotation',hAngle);
end

function reference = full_matrix_difference(q,Q,setup,opt,h)
%FULL_MATRIX_DIFFERENCE checks direct columns against differentiated R_SU*Mbody.
P = size(q,1);
Dbody = zeros(3,3,P);
E = zeros(3,3,5);
E(:,:,1) = diag([1,-1,0])/sqrt(2);
E(:,:,2) = diag([1,1,-2])/sqrt(6);
E(1,2,3)=1/sqrt(2); E(2,1,3)=E(1,2,3);
E(1,3,4)=1/sqrt(2); E(3,1,4)=E(1,3,4);
E(2,3,5)=1/sqrt(2); E(3,2,5)=E(2,3,5);
timer = tic;
maxSurfaceResidual = 0;
for j = 1:6*P
    coordinate = mod(j-1,6)+1;
    if coordinate>3, continue; end
    body = floor((j-1)/6)+1;
    columns = cell(2,1);
    for side = 1:2
        shiftedQ = q;
        shiftedQ(body,coordinate) = shiftedQ(body,coordinate)+(2*side-3)*h;
        maps = measure_maps(shiftedQ,Q,setup,opt);
        maxSurfaceResidual = max(maxSurfaceResidual,maps.maxSurfaceResidual);
        columns{side} = maps.A(:,j);
    end
    columnDerivative = (columns{2}-columns{1})/(2*h);
    for p = 1:P
        for alpha = 1:5
            Dbody(:,:,p) = Dbody(:,:,p)+ ...
                columnDerivative(5*(p-1)+alpha)*E(:,:,alpha);
        end
    end
end
reference = struct('Dbody',Dbody, ...
    'resistanceRightHandSides',12*6*P, ...
    'mobilityRightHandSides',12*6*P, ...
    'maxSurfaceResidual',maxSurfaceResidual,'seconds',toc(timer));
end
