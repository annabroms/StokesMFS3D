%% OCT23_STAGNANT_CAP: full stagnant-cap solver and convergence experiment
% Investigates how boundary residuals, drag, and settling speed change as
% the source and collocation resolutions grow at one fixed cap angle.
% The configured sequence also changes ringLevels, collocation ratios, and
% optionally ring azimuth counts; read those arrays together when comparing
% resolutions. ringLevels counts added collocation rings in each region.
% The cap is theta<phi from +z; gravity is axial.
% At each resolution, solve both prescribed-speed and prescribed-force
% problems, then check note conditions (30)-(33) at independent surface
% points and compare drag and speed with equation (50).
% Also show each collocation geometry, surface residual maps, the finest-geometry
% geometry and drop-relative velocity, and optional limits (51)-(53).
% Edit the controls below. Generated figures stay open; PNG and CSV copies
% are written to experiments/results/.

clear; close all; 

experimentDirectory = fileparts(mfilename('fullpath'));
repositoryDirectory = fileparts(fileparts(experimentDirectory));
addpath(fullfile(repositoryDirectory,'drops'), ...
    fullfile(repositoryDirectory,'images'), ...
    fullfile(repositoryDirectory,'common','dense'), ...
    fullfile(repositoryDirectory,'geometry','spheres'), ...
    fullfile(repositoryDirectory,'tools'));
outputDirectory = fullfile(experimentDirectory,'results');
if ~isfolder(outputDirectory), mkdir(outputDirectory); end

% Edit this configuration to change the physical or numerical problem.
useExtraCollocation = true; % Set false for base spherical-design nodes only.
options = stagnant_cap_options(struct( ...
    'radius',1,'muOut',1,'muIn',2,'phi',pi/4, ...
    'innerProxyRatio',0.63,'outerProxyRatio',2, ...
    'svdTolerance',1e-12,'collocationRatio',1.5, ...
    'useExtraCollocation',useExtraCollocation, ...
    'gravity',[0 0 -1],'densityDifference',1));
% For a smaller trial, set sourceRequests to [50 100 200].
sourceRequests = [500 700 1000];
ringLevels = [3 4 5];
collocationRatios = [1.5 1.5 2.0];
ringAzimuthCounts = [NaN NaN 25];
prescribedSpeed = -1;
runLimitCases = true;
assert(numel(sourceRequests)==numel(ringLevels));
Fg = (4*pi*options.radius^3/3)* ...
    options.densityDifference*options.gravity;
assert(Fg(1)==0 && Fg(2)==0 && Fg(3)~=0);
exactDrag = stagnant_cap_exact_drag( ...
    options.phi,options.muIn,options.muOut,options.radius);
exactSpeed = Fg(3)/exactDrag;

grids = cell(numel(sourceRequests),1);
caseOptions = cell(numel(sourceRequests),1);
for k = 1:numel(sourceRequests)
    o = options;
    o.sourceRequest = sourceRequests(k);
    o.ringLevels = ringLevels(k);
    o.collocationRatio = collocationRatios(k);
    if isfinite(ringAzimuthCounts(k))
        o.ringAzimuthCount = ringAzimuthCounts(k);
    end
    caseOptions{k} = stagnant_cap_options(o);
    grids{k} = stagnant_cap_geometry(caseOptions{k});
end
maxM = max(cellfun(@(g) g.collocationCount,grids));
% Keep near-edge validation points even when extra collocation is disabled.
finestOffset = options.edgeBandAngle* ...
    2^(-max(1,max(ringLevels)-1));
validation = stagnant_cap_validation_grid(options,maxM,finestOffset);
fprintf('Shared independent validation: %d points, minimum edge distance %.6g.\n', ...
    validation.count,validation.minimumEdgeDistance);

requestedN = zeros(numel(grids),1);
actualN = requestedN; collocationM = requestedN;
rankSpeed = requestedN; rankForce = requestedN;
svdRows = requestedN; svdCols = requestedN;
svdInputMB = requestedN; svdKnownDenseMB = requestedN;
drag = requestedN; speed = requestedN;
dragError = requestedN; speedError = requestedN;
speedFarMax = requestedN; speedNearMax = requestedN;
forceFarMax = requestedN; forceNearMax = requestedN;
forceResidual = requestedN; transverseDrift = requestedN;
torqueMagnitude = requestedN;
allValidation = table();
for k = 1:numel(grids)
    o = caseOptions{k};
    geometry = grids{k};
    collocationFigure = plot_stagnant_cap_collocation(geometry,o,'on');
    collocationName = sprintf('oct23_stagnant_cap_collocation_N%04d.png', ...
        geometry.N);
    exportgraphics(collocationFigure, ...
        fullfile(outputDirectory,collocationName),'Resolution',180);
    prescribed = solve_stagnant_cap(o, ...
        struct('kind','speed','Uz',prescribedSpeed),geometry);
    driven = solve_stagnant_cap(o,struct('kind','force','Fg',Fg),geometry);
    checkSpeed = stagnant_cap_validate(prescribed,validation,o);
    checkForce = stagnant_cap_validate(driven,validation,o);
    % Use the same independent points; each condition scales its colourbar.
    maps = {{prescribed,checkSpeed,'speed'}, ...
        {driven,checkForce,'force'}};
    for mapIndex = 1:numel(maps)
        entry = maps{mapIndex};
        mapFigure = plot_stagnant_cap_residuals( ...
            entry{2},validation,entry{1},o,'on');
        mapName = sprintf('oct23_stagnant_cap_residuals_%s_N%04d.png', ...
            entry{3},geometry.N);
        exportgraphics(mapFigure,fullfile(outputDirectory,mapName), ...
            'Resolution',150);
    end
    fprintf('\nN requested/actual %d/%d, M=%d, closest collocation distance %.5g rad.\n', ...
        sourceRequests(k),geometry.N,geometry.collocationCount, ...
        min(geometry.edgeDistance)/o.radius);
    fprintf('Prescribed speed: U_z=%.6g, drag=%.8g, exact=%.8g.\n', ...
        prescribedSpeed,prescribed.sourceForce(3)/prescribedSpeed,exactDrag);
    fprintf('Prescribed force: F_g,z=%.6g, U_z=%.8g, exact=%.8g.\n', ...
        Fg(3),driven.U(3),exactSpeed);
    fprintf('SVD ranks speed/force %d/%d of %d; cutoff %.1e.\n', ...
        prescribed.retainedRank,driven.retainedRank, ...
        6*geometry.N,o.svdTolerance);
    speedRows = checkSpeed.table;
    forceRows = checkForce.table;
    speedRows.caseIndex = k*ones(height(speedRows),1);
    speedRows.mode = repmat("speed",height(speedRows),1);
    forceRows.caseIndex = k*ones(height(forceRows),1);
    forceRows.mode = repmat("force",height(forceRows),1);
    caseValidation = [speedRows;forceRows];
    fprintf('Independent validation residuals for N=%d:\n',geometry.N);
    disp(caseValidation(:,{'mode','condition','region','zone', ...
        'maximum','rms'}));
    allValidation = [allValidation;caseValidation]; %#ok<AGROW>
    requestedN(k) = sourceRequests(k);
    actualN(k) = geometry.N;
    collocationM(k) = geometry.collocationCount;
    rankSpeed(k) = prescribed.retainedRank;
    rankForce(k) = driven.retainedRank;
    svdRows(k) = prescribed.matrixSize(1);
    svdCols(k) = prescribed.matrixSize(2);
    svdInputMB(k) = prescribed.svdInputMB;
    svdKnownDenseMB(k) = prescribed.svdKnownDenseMB;
    assert(isequal(prescribed.matrixSize,driven.matrixSize));
    drag(k) = prescribed.sourceForce(3)/prescribedSpeed;
    speed(k) = driven.U(3);
    dragError(k) = abs(drag(k)/exactDrag-1);
    speedError(k) = abs(speed(k)/exactSpeed-1);
    speedFarMax(k) = max(speedRows.maximum(speedRows.zone=="elsewhere"));
    speedNearMax(k) = max(speedRows.maximum(speedRows.zone=="near edge"));
    forceFarMax(k) = max(forceRows.maximum(forceRows.zone=="elsewhere"));
    forceNearMax(k) = max(forceRows.maximum(forceRows.zone=="near edge"));
    forceResidual(k) = driven.forceResidual/norm(Fg);
    transverseDrift(k) = norm(driven.U(1:2))/abs(driven.U(3));
    torqueMagnitude(k) = norm(driven.torque)/(norm(Fg)*o.radius);
    if k==numel(grids)
        finestSolution = driven;
        finestGrid = geometry;
        finestOptions = o;
    end
end
summary = table(requestedN,actualN,collocationM,svdRows,svdCols, ...
    svdInputMB,svdKnownDenseMB,rankSpeed,rankForce, ...
    drag,dragError,speed,speedError,speedFarMax,speedNearMax, ...
    forceFarMax,forceNearMax,forceResidual,transverseDrift,torqueMagnitude);
summary.ringLevels = ringLevels(:);
disp(summary);
writetable(summary,fullfile(outputDirectory,'oct23_stagnant_cap_summary.csv'));
writetable(allValidation,fullfile(outputDirectory, ...
    'oct23_stagnant_cap_validation.csv'));

[geometryFigure,closeFigure] = ...
    plot_stagnant_cap_geometry(finestGrid,finestOptions,'on');
exportgraphics(geometryFigure,fullfile(outputDirectory, ...
    'oct23_stagnant_cap_geometry_3d.png'),'Resolution',180);
exportgraphics(closeFigure,fullfile(outputDirectory, ...
    'oct23_stagnant_cap_geometry_edge.png'),'Resolution',180);
[surfaceFigure,profileFigure,edgeDiagnostic] = ...
    plot_stagnant_cap_velocity(finestSolution,finestOptions,'on');
exportgraphics(surfaceFigure,fullfile(outputDirectory, ...
    'oct23_stagnant_cap_surface_velocity.png'),'Resolution',180);
exportgraphics(profileFigure,fullfile(outputDirectory, ...
    'oct23_stagnant_cap_edge_profile.png'),'Resolution',180);
fprintf('Square-root fit coefficient %.6g over s in [%.6g, %.6g].\n', ...
    edgeDiagnostic.fitC,edgeDiagnostic.fitWindow);
fprintf('Closest plotted mobile value: s=%.6g, relative u_theta=%.6g; theoretical edge limit is zero.\n', ...
    edgeDiagnostic.closestMobileS,edgeDiagnostic.closestMobileVelocity);

if runLimitCases
    % Equation (51): clean drop; (52): immobilized sphere;
    % (53): bubble at a finite cap angle, approached with muIn>0.
    labels = ["clean drop";"immobilized sphere";"bubble limit"];
    angles = [0;pi;pi/2];
    interiorViscosities = [options.muIn;options.muIn;1e-2*options.muOut];
    numericalDrag = zeros(3,1); referenceDrag = zeros(3,1);
    numericalSpeed = zeros(3,1); referenceSpeed = zeros(3,1);
    limitSvdRows = zeros(3,1); limitSvdCols = zeros(3,1);
    limitSvdInputMB = zeros(3,1); limitSvdKnownDenseMB = zeros(3,1);
    for j = 1:3
        o = options;
        o.phi = angles(j);
        o.muIn = interiorViscosities(j);
        o.sourceRequest = 100;
        o.ringLevels = 4;
        o.edgeBandAngle = [];
        o = stagnant_cap_options(o);
        g = stagnant_cap_geometry(o);
        r = solve_stagnant_cap(o, ...
            struct('kind','speed','Uz',prescribedSpeed),g);
        m = solve_stagnant_cap(o,struct('kind','force','Fg',Fg),g);
        numericalDrag(j) = r.sourceForce(3)/prescribedSpeed;
        referenceDrag(j) = stagnant_cap_exact_drag( ...
            o.phi,o.muIn,o.muOut,o.radius);
        numericalSpeed(j) = m.U(3);
        referenceSpeed(j) = Fg(3)/referenceDrag(j);
        limitSvdRows(j) = r.matrixSize(1);
        limitSvdCols(j) = r.matrixSize(2);
        limitSvdInputMB(j) = r.svdInputMB;
        limitSvdKnownDenseMB(j) = r.svdKnownDenseMB;
        assert(isequal(r.matrixSize,m.matrixSize));
    end
    relativeDragError = abs(numericalDrag./referenceDrag-1);
    relativeSpeedError = abs(numericalSpeed./referenceSpeed-1);
    limits = table(labels,angles,interiorViscosities, ...
        limitSvdRows,limitSvdCols,limitSvdInputMB, ...
        limitSvdKnownDenseMB,numericalDrag, ...
        referenceDrag,relativeDragError,numericalSpeed,referenceSpeed, ...
        relativeSpeedError);
    disp(limits);
    writetable(limits,fullfile(outputDirectory, ...
        'oct23_stagnant_cap_limits.csv'));
    bubbleExact = stagnant_cap_exact_drag(pi/2,0, ...
        options.muOut,options.radius);
    fprintf('Equation (53) zero-viscosity bubble reference: %.8g; finite-ratio reference %.8g.\n', ...
        bubbleExact,referenceDrag(3));
end
