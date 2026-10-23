% OCT15_2026_TWO_DROP_AXIAL_PRESSURE_GAP_SWEEP
% Verify exterior pressure with/without images at delta=0.1, then sweep
% two proxy resolutions with images enabled for near pairs. Pressure uses
% p_0 -> 0 at infinity; interior pressure gauges never enter these plots.
% Write one final MAT file after the sweep and pressure-field evaluation.
% Replot later with plot_drop_pressure_sweep(outputFile).
clear; close all

%% Physical problem and the two proxy resolutions
radius = 1;
innerProxyRadius = 0.63;
outerProxyRadius = 2;
S1 = 2;
interiorViscosity = 10;
e = repmat([0 0 1],2,1);          % Active directions along the line of centres. z direction gives much larger errors than x direction
e = repmat([1 0 0],2,1);
gapList = [1 0.5 0.2 0.1 0.05 0.02 0.01 0.005 0.002 0.001];
gapList = [0.1];
representativeGap = 0.1;
requestedN = [700 1000];            % Requested nodes per proxy, per drop.
collocationToSourceRatio = 2;

%% Image and cap controls, using the two-drop demo's source mixture
N_image = 30;                     % [] selects the shared adaptive gap rule.
% imageFamilies tokens: S Stokeslet, R rotlet, D potential dipole, T legacy
% stresslet, T_t line-parallel, T_n first transverse, T_n2 second transverse
% stresslet, Q_s trace-free axial Stokeslet quadrupole. T excludes T_*.
imageFamilies = "SDRT_nT_t";
mirror_ext = true;
nearThreshold = 0.15;             % Images only for delta/radius below this value.
outerImageOffsetFraction = 0.001;
capMode = "exponential";          % Or two overlapping "fibonacci" caps.
capAngles = [pi/5 pi/60];
capPointCounts = [800 200];        % Exponential mode combines the two counts.
capExponent = 3;
capLayers = 18;
useRowScaling = true;
plot_grid = false;                % If true, show one drop before every solve.

%% Accuracy controls and independent validation
svdCutoff = 1e-12;
gmresTolerance = 1e-10;
maxIterations = 350;
validationResolution = 4000;
validationGapCounts = [4001 1001]; % Wide and narrow independent Fibonacci caps.
checkGap = 0.1;
pressureCheckTolerance = 1e-3;     % Relative infinity-norm agreement, not equality.
pressureReferenceN = 1000;         % Additional no-image refinement at checkGap.
batchSize = 128;
planeGridResolution = [401 401];   % [x,z] nodes in the zoomed pressure window.
gapZoomFactor = 3;                % Vertical half-width in sqrt(radius*delta).
view_solution = true;

%% Paths, physical scales, and common solver options
experimentDirectory = fileparts(mfilename('fullpath'));
repositoryDirectory = fileparts(fileparts(experimentDirectory));
addpath(genpath(repositoryDirectory));
outputDirectory = fullfile(experimentDirectory,'results','oct19_2026_axial_pressure_images');
gapList = sort(gapList,'descend');
if numel(requestedN)~=2 || ~ismember(representativeGap,gapList)
    error('Choose two resolutions and include representativeGap in gapList.');
end
a = radius*ones(2,1);
S = S1*ones(2,1);
e = e./vecnorm(e,2,2);
velocityScale = max(abs(S1)*radius,realmin);
tractionScale = max(abs(S1),realmin);
residualScale = [velocityScale velocityScale tractionScale];
base = struct('N',max(requestedN),'proxyRadius',innerProxyRadius/radius, ...
    'outerProxyRatio',outerProxyRadius/radius,'collocationRatio',collocationToSourceRatio, ...
    'useImages',true,'N_image',N_image,'imageFamilies',imageFamilies, ... % See token legend above.
    'mirror_ext',mirror_ext,'nearThreshold',nearThreshold, ...
    'outerImageOffsetFraction',outerImageOffsetFraction, ...
    'capMode',capMode,'capAngles',capAngles,'capPointCounts',capPointCounts, ...
    'capExponent',[capExponent capExponent],'capLayers',[capLayers capLayers], ...
    'useRowScaling',useRowScaling,'rowWeightMode',"local", ...
    'interiorViscosity',interiorViscosity,'completionScale',1/(8*pi*radius), ...
    'velocityScale',velocityScale,'tractionScale',tractionScale, ...
    'tractionRowWeight',1,'scale_columns',false,'svdCutoff',svdCutoff, ...
    'gmresTolerance',gmresTolerance,'maxIterations',maxIterations, ...
    'batchSize',batchSize,'solverMode',"oneBody",'plot_grid',plot_grid);
base = drop_image_options(base,2,a);

%% First verify the pressure field at delta=0.1 on identical exterior targets
q = [-(radius+checkGap/2) 0 0;radius+checkGap/2 0 0];
checkGrid = drop_image_discretization(q,a,base);
[checkPoints,checkNormals] = drop_image_validation_grid( ...
    checkGrid,validationResolution,validationGapCounts,capAngles);
% Check the whole exterior cross-section, a finer gap patch, and both axes.
[x,z] = meshgrid(linspace(-3*radius,3*radius,61),linspace(-2*radius,2*radius,61));
[xg,zg] = meshgrid(linspace(-2*checkGap,2*checkGap,61), ...
    linspace(-sqrt(radius*checkGap),sqrt(radius*checkGap),81));
axialX = linspace(-0.49*checkGap,0.49*checkGap,101).';
transverseZ = linspace(-radius,radius,201).';
checkTargets = [x(:),zeros(numel(x),1),z(:);xg(:),zeros(numel(xg),1),zg(:); ...
    axialX,zeros(numel(axialX),2);zeros(numel(transverseZ),2),transverseZ];
outside = true(size(checkTargets,1),1);
for i = 1:2
    outside = outside & sum((checkTargets-q(i,:)).^2,2)>radius^2;
end
checkTargets = checkTargets(outside,:);
checkN = [base.N;base.N;pressureReferenceN];
checkImages = [false;true;false];
checkActualN = zeros(3,1);
checkPressure = zeros(size(checkTargets,1),3);
checkMidpoint = zeros(3,1);
checkResiduals = zeros(3,3);
checkIterations = zeros(3,1);
checkOriginal = zeros(3,1);
for k = 1:3
    opt = base; opt.useImages = checkImages(k); opt.N = checkN(k);
    fprintf('\nPressure verification at delta=%.3g, N=%d, useImages=%d\n', ...
        checkGap,opt.N,opt.useImages);
    [U,lambda,gamma,g,info] = solve_drops_with_images(q,a,S,e,opt);
    % grid.body.basis contains every exterior-flow family, including both T's.
    % Outer sources affect this pressure through the solve, not by direct summation.
    checkPressure(:,k) = image_pressure(g,lambda,checkTargets,1);
    checkMidpoint(k) = image_pressure(g,lambda,[0 0 0],1);
    maxima = drop_image_validate_surface_conditions( ...
        g,lambda,gamma,U,S,e,info.options,checkPoints,checkNormals);
    checkResiduals(k,:) = max(maxima,[],1)./residualScale;
    checkIterations(k) = info.iterations;
    checkOriginal(k) = info.originalResidual;
    checkActualN(k) = g.N;
    if k==2
        % Reuse the pressure-check image solution when the sweep reaches this case.
        cachedGrid = g; cachedLambda = lambda; cachedGamma = gamma;
        cachedU = U; cachedInfo = info; cachedMaxima = maxima;
    end
end
pressureDifference = norm(checkPressure(:,2)-checkPressure(:,1),inf);
pressureRelative = pressureDifference/max(norm(checkPressure(:,1),inf),realmin);
midpointRelative = abs(checkMidpoint(2)-checkMidpoint(1))/max(abs(checkMidpoint(1)),realmin);
refinedRelative = norm(checkPressure(:,3)-checkPressure(:,2),inf)/ ...
    max(norm(checkPressure(:,2),inf),realmin);
pressureCheckPassed = pressureRelative<=pressureCheckTolerance && ...
    midpointRelative<=pressureCheckTolerance;
checkTable = table(checkN,checkActualN,checkImages,checkMidpoint,checkIterations,checkOriginal, ...
    checkResiduals(:,1),checkResiduals(:,2),checkResiduals(:,3), ...
    'VariableNames',{'NRequested','NActual','Images','MidpointPressure','GMRESIterations', ...
    'OriginalResidual','ContinuityRel','NormalRel','TractionRel'});
disp(checkTable);
fprintf('Pressure difference: absolute max %.3e; relative field %.3e; relative midpoint %.3e; tolerance %.3e\n', ...
    pressureDifference,pressureRelative,midpointRelative,pressureCheckTolerance);
fprintf('Refined no-image field versus image field: relative difference %.3e.\n',refinedRelative);
if ~pressureCheckPassed
    warning('axial_pressure_gap_sweep:PressureCheck', ...
        ['The no-image and image pressures do not agree at the requested tolerance. ', ...
        'Retain this failed check when interpreting the image sweep.']);
else
    fprintf('Pressure verification passed at the requested tolerance.\n');
end
fprintf('Starting the two-resolution image sweep; independent residuals are reported for every case.\n');

% A transverse gap profile makes any pressure disagreement visible.
if view_solution
    rows = find(abs(checkTargets(:,1))<eps & abs(checkTargets(:,3))<=radius);
    [~,order] = sort(checkTargets(rows,3)); rows = rows(order);
    figure('Name','Pressure verification at delta=0.1','Color','w');
    hold on
    labels = ["without images","with images","without images"];
    for k = 1:3
        plot(checkTargets(rows,3),checkPressure(rows,k),'LineWidth',1.3, ...
            'DisplayName',sprintf('N=%d, %s',checkActualN(k),labels(k)));
    end
    xlabel('z at x=y=0'); ylabel('Exterior pressure p_0');
    title(sprintf('Pressure comparison at delta=%.3g',checkGap));
    grid on; legend('Location','best');
end

%% Sweep gaps with a common independent grid at each separation
names = {'NRequested','NActual','Gap','ImagesPerLine','CollocationPerDrop', ...
    'U1Toward','U2Toward','GapRate','MidpointPressure','PressureMagnitude', ...
    'MaxTransverseSpeed','GMRESIterations','GMRESResidual','OriginalResidual', ...
    'ContinuityRel','NormalRel','TractionRel','Rank1','Rank2', ...
    'SetupSeconds','SolveSeconds','ForceSumNorm'};
sweep = array2table(nan(2*numel(gapList),numel(names)),'VariableNames',names);
representativeGrids = cell(2,1);
representativeStrengths = cell(2,1);
row = 0;
for resolution = 1:2
    opt = base; opt.N = requestedN(resolution);
    for delta = gapList
        row = row+1;
        q = [-(radius+delta/2) 0 0;radius+delta/2 0 0];
        fprintf('\nSweep %d/%d: requested N=%d, delta=%.4g\n', ...
            row,height(sweep),opt.N,delta);
        if opt.N==base.N && delta==checkGap
            U = cachedU; lambda = cachedLambda; gamma = cachedGamma;
            g = cachedGrid; info = cachedInfo; maxima = cachedMaxima;
            fprintf('Reusing the pressure-check image solution at delta=0.1.\n');
        else
            [U,lambda,gamma,g,info] = solve_drops_with_images(q,a,S,e,opt);
            [points,normals] = drop_image_validation_grid( ...
                g,validationResolution,validationGapCounts,capAngles);
            maxima = drop_image_validate_surface_conditions( ...
                g,lambda,gamma,U,S,e,info.options,points,normals);
        end
        pMid = image_pressure(g,lambda,[0 0 0],1);
        boundaryMaxima = max(maxima,[],1)./residualScale;
        sweep{row,:} = [opt.N,g.N,delta,g.body(1).N_image, ...
            size(g.body(1).collocation,1),U(1,1),-U(2,1),U(2,1)-U(1,1), ...
            pMid,abs(pMid),max(vecnorm(U(:,2:3),2,2)),info.iterations, ...
            info.regularizedResidual,info.originalResidual,boundaryMaxima, ...
            info.retainedRanks(:).',info.setupTime,info.solveTime,max(info.forceSumNorm)];
        fprintf('U1x=% .5e, U2x=% .5e, gap rate=% .5e, p_mid=% .5e; GMRES %d, BC [%.2e %.2e %.2e]\n', ...
            U(1,1),U(2,1),U(2,1)-U(1,1),pMid,info.iterations,boundaryMaxima);
        if ~info.converged
            warning('axial_pressure_gap_sweep:GMRES', ...
                'GMRES missed its tolerance at N=%d, delta=%.4g.',opt.N,delta);
        end
        if delta==representativeGap
            representativeGrids{resolution} = g;
            representativeStrengths{resolution} = lambda;
        end
    end
end
sweep.GMRESConverged = sweep.GMRESResidual<=gmresTolerance;
sweep.IterationLimitReached = sweep.GMRESIterations>=maxIterations;
fprintf('\nExterior pressure has p_0 -> 0 at infinity. Images activate for delta/a < %.3g.\n',nearThreshold);
disp(sweep);
% Resolution differences are diagnostics, not an extrapolated error estimate.
coarse = sweep(sweep.NRequested==requestedN(1),:);
fine = sweep(sweep.NRequested==requestedN(2),:);
pressureChange = abs(fine.MidpointPressure-coarse.MidpointPressure)./ ...
    max(abs(fine.MidpointPressure),realmin);
resolutionTable = table(fine.Gap,coarse.MidpointPressure,fine.MidpointPressure, ...
    pressureChange,'VariableNames',{'Gap','PressureCoarse','PressureFine','RelativeChange'});
disp(resolutionTable);

%% Zoomed exterior pressure: common colour limits and exactly white at zero
zExtent = min(0.8*radius,gapZoomFactor*sqrt(radius*representativeGap));
xExtent = 1.25*(radius+representativeGap/2-sqrt(radius^2-zExtent^2));
xGrid = linspace(-xExtent,xExtent,planeGridResolution(1));
zGrid = linspace(-zExtent,zExtent,planeGridResolution(2));
[gridX,gridZ] = meshgrid(xGrid,zGrid);
planeTargets = [gridX(:),zeros(numel(gridX),1),gridZ(:)];
q = representativeGrids{1}.centres;
exterior = true(size(planeTargets,1),1);
for i = 1:2
    exterior = exterior & sum((planeTargets-q(i,:)).^2,2)>radius^2;
end
planePressure = nan(size(planeTargets,1),2);
for resolution = 1:2
    planePressure(exterior,resolution) = image_pressure( ...
        representativeGrids{resolution},representativeStrengths{resolution}, ...
        planeTargets(exterior,:),1);
end
pressureLimit = max(max(abs(planePressure(exterior,:)),[],'all'),eps);
pressureColours = interp1([0 0.5 1], ...
    [0.10 0.27 0.60;1 1 1;0.65 0.07 0.12],linspace(0,1,257).');

%% Save the completed sweep, checks, and representative pressure field once
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
outputFile = fullfile(outputDirectory,'gap_sweep.mat');
save(outputFile,'sweep','resolutionTable','requestedN','gapList','base', ...
    'radius','S','e','validationResolution','validationGapCounts', ...
    'checkTable','checkGap','pressureDifference','pressureRelative', ...
    'midpointRelative','refinedRelative','pressureCheckPassed','pressureCheckTolerance', ...
    'representativeGap','representativeGrids','representativeStrengths', ...
    'xGrid','zGrid','planePressure','exterior','pressureLimit','pressureColours','q');
fprintf('Final results saved in %s\n',outputFile);

% Use the same plotter for live results and later inspection of the saved file.
if view_solution
    plot_drop_pressure_sweep(outputFile);
end
