%% OCT23_STAGNANT_CAP_N_SWEEP: accuracy versus source count at fixed cap angle
% Investigates how the four imposed boundary-condition maxima and the
% equation (50) drag/speed errors change as N grows on each uniform proxy.
% Hold phi, physical parameters, and proxy radii fixed. The ringLevels
% array below counts added collocation rings per region and changes with N
% by default, so this is a coupled source
% and edge-collocation refinement sequence; set it constant to isolate N.
% Use the same independent validation points for every N and both solve
% modes. The six-panel figure shows mobile continuity, mobile no-penetration,
% mobile traction jump, cap tangential velocity, drag error, and speed error.
% Edit sourceRequests, ringLevels, phi, and useExtraCollocation below.
% Figure and summary/validation CSV files go to experiments/results/.

useExtraCollocation = true; % Set false for base spherical-design nodes only.
phi = pi/2;
sourceRequests = [50 100 200 500 600 700 800 1000]; % Minimum N on EACH proxy surface.
ringLevels = [3 4 5 5 5 5 5 5]; % Used for graded collocation and validation spacing.
prescribedSpeed = -1;
experimentDirectory = fileparts(mfilename('fullpath'));
repositoryDirectory = fileparts(fileparts(experimentDirectory));
addpath(fullfile(repositoryDirectory,'drops'), ...
    fullfile(repositoryDirectory,'images'), ...
    fullfile(repositoryDirectory,'common','dense'), ...
    fullfile(repositoryDirectory,'geometry','spheres'), ...
    fullfile(repositoryDirectory,'tools'));
outputDirectory = fullfile(experimentDirectory,'results');
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
options = stagnant_cap_options(struct( ...
    'radius',1,'muOut',1,'muIn',2,'phi',phi, ...
    'innerProxyRatio',0.4,'outerProxyRatio',2, ...
    'collocationRatio',1.5,'svdTolerance',1e-12, ...
    'useExtraCollocation',useExtraCollocation, ...
    'gravity',[0 0 -1],'densityDifference',1));
assert(options.phi>0 && options.phi<pi, ...
    'The four-region diagnostic requires 0 < phi < pi.');
assert(numel(sourceRequests)==numel(ringLevels) && ...
    all(diff(sourceRequests)>0));
assert(prescribedSpeed~=0 && isfinite(prescribedSpeed));
Fg = (4*pi*options.radius^3/3)* ...
    options.densityDifference*options.gravity;
exactDrag = stagnant_cap_exact_drag(options.phi, ...
    options.muIn,options.muOut,options.radius);
exactSpeed = Fg(3)/exactDrag;
assert(exactSpeed~=0);

caseCount = numel(sourceRequests);
caseOptions = cell(caseCount,1);
geometries = cell(caseCount,1);
for k = 1:caseCount
    o = options;
    o.sourceRequest = sourceRequests(k);
    o.ringLevels = ringLevels(k);
    caseOptions{k} = stagnant_cap_options(o);
    geometries{k} = stagnant_cap_geometry(caseOptions{k});
end
maxM = max(cellfun(@(g) g.collocationCount,geometries));
% Validation stays edge-refined in both collocation modes and is shared
% across N. Its finest ring is half this nominal finest collocation offset.
finestOffset = options.edgeBandAngle* ...
    2^(-max(1,max(ringLevels)-1));
validation = stagnant_cap_validation_grid(options,maxM,finestOffset);
assert(any(validation.cap) && any(validation.mobile));
fprintf(['Independent validation: %d points (>=5 x largest M), ', ...
    'minimum edge distance %.6g a.\n'], ...
    validation.count,validation.minimumEdgeDistance/options.radius);
fprintf(['Velocity residuals use U_*; traction residuals use ', ...
    'mu_out U_*/a. Maxima include near-edge validation points.\n']);

requestedN = sourceRequests(:);
actualN = zeros(caseCount,1);
collocationM = zeros(caseCount,1);
edgeRingNodeCount = zeros(caseCount,1);
collocationMinimumEdgeDistance = zeros(caseCount,1);
svdRows = zeros(caseCount,1); svdCols = zeros(caseCount,1);
svdInputMB = zeros(caseCount,1);
svdKnownDenseMB = zeros(caseCount,1);
rankSpeed = zeros(caseCount,1); rankForce = zeros(caseCount,1);
velocityScaleSpeed = zeros(caseCount,1);
velocityScaleForce = zeros(caseCount,1);
tractionScaleSpeed = zeros(caseCount,1);
tractionScaleForce = zeros(caseCount,1);
drag = zeros(caseCount,1); speed = zeros(caseCount,1);
dragError = zeros(caseCount,1); speedError = zeros(caseCount,1);
% Columns: mobile continuity, mobile normal velocity, mobile traction,
% cap exterior tangential velocity. Both solve modes are reported.
residualSpeed = zeros(caseCount,4);
residualForce = zeros(caseCount,4);
allValidation = table();

for k = 1:caseCount
    o = caseOptions{k};
    geometry = geometries{k};
    resistance = solve_stagnant_cap(o, ...
        struct('kind','speed','Uz',prescribedSpeed),geometry);
    mobility = solve_stagnant_cap(o, ...
        struct('kind','force','Fg',Fg),geometry);
    checkSpeed = stagnant_cap_validate(resistance,validation,o);
    checkForce = stagnant_cap_validate(mobility,validation,o);
    residualSpeed(k,:) = stagnant_cap_boundary_maxima( ...
        checkSpeed,validation,resistance);
    residualForce(k,:) = stagnant_cap_boundary_maxima( ...
        checkForce,validation,mobility);
    actualN(k) = geometry.N;
    collocationM(k) = geometry.collocationCount;
    edgeRingNodeCount(k) = geometry.collocationCount-geometry.baseCount;
    collocationMinimumEdgeDistance(k) = min(geometry.edgeDistance);
    svdRows(k) = resistance.matrixSize(1);
    svdCols(k) = resistance.matrixSize(2);
    svdInputMB(k) = resistance.svdInputMB;
    svdKnownDenseMB(k) = resistance.svdKnownDenseMB;
    assert(isequal(resistance.matrixSize,mobility.matrixSize));
    rankSpeed(k) = resistance.retainedRank;
    rankForce(k) = mobility.retainedRank;
    velocityScaleSpeed(k) = resistance.velocityScale;
    velocityScaleForce(k) = mobility.velocityScale;
    tractionScaleSpeed(k) = resistance.tractionScale;
    tractionScaleForce(k) = mobility.tractionScale;
    drag(k) = -resistance.FH(3)/prescribedSpeed;
    speed(k) = mobility.U(3);
    dragError(k) = abs(drag(k)/exactDrag-1);
    speedError(k) = abs(speed(k)/exactSpeed-1);
    fprintf('\nN=%d, M=%d, extra collocation=%d:\n', ...
        actualN(k),collocationM(k),useExtraCollocation);
    fprintf('  mobile max [continuity, normal, traction], cap max [tangent]\n');
    fprintf('  speed solve: %.3g  %.3g  %.3g  %.3g\n',residualSpeed(k,:));
    fprintf('  force solve: %.3g  %.3g  %.3g  %.3g\n',residualForce(k,:));
    fprintf('  drag error %.3g, speed error %.3g\n', ...
        dragError(k),speedError(k));
    speedRows = checkSpeed.table;
    forceRows = checkForce.table;
    speedRows.actualN = repmat(actualN(k),height(speedRows),1);
    speedRows.mode = repmat("speed",height(speedRows),1);
    forceRows.actualN = repmat(actualN(k),height(forceRows),1);
    forceRows.mode = repmat("force",height(forceRows),1);
    allValidation = [allValidation;speedRows;forceRows]; %#ok<AGROW>
end

extraCollocation = options.useExtraCollocation & true(caseCount,1);
validationCount = repmat(validation.count,caseCount,1);
validationMinimumEdgeDistance = repmat( ...
    validation.minimumEdgeDistance,caseCount,1);
referenceDrag = repmat(exactDrag,caseCount,1);
referenceSpeed = repmat(exactSpeed,caseCount,1);
summary = table(requestedN,actualN,collocationM,extraCollocation, ...
    edgeRingNodeCount,collocationMinimumEdgeDistance,validationCount, ...
    validationMinimumEdgeDistance,svdRows,svdCols,svdInputMB, ...
    svdKnownDenseMB,rankSpeed,rankForce,velocityScaleSpeed, ...
    velocityScaleForce,tractionScaleSpeed,tractionScaleForce, ...
    drag,referenceDrag,dragError,speed,referenceSpeed,speedError);
summary.ringLevels = ringLevels(:);
summary.mobileContinuitySpeed = residualSpeed(:,1);
summary.mobileContinuityForce = residualForce(:,1);
summary.mobileNormalSpeed = residualSpeed(:,2);
summary.mobileNormalForce = residualForce(:,2);
summary.mobileTractionSpeed = residualSpeed(:,3);
summary.mobileTractionForce = residualForce(:,3);
summary.capTangentialSpeed = residualSpeed(:,4);
summary.capTangentialForce = residualForce(:,4);
disp(summary(:,{'actualN','collocationM', ...
    'mobileContinuitySpeed','mobileNormalSpeed', ...
    'mobileTractionSpeed','capTangentialSpeed', ...
    'mobileContinuityForce','mobileNormalForce', ...
    'mobileTractionForce','capTangentialForce', ...
    'dragError','speedError'}));

if options.useExtraCollocation, modeName = 'graded'; else, modeName = 'base'; end
stem = sprintf('oct23_stagnant_cap_N_sweep_phi%05.1f_%s', ...
    rad2deg(options.phi),modeName);
writetable(summary,fullfile(outputDirectory,[stem '.csv']));
writetable(allValidation,fullfile(outputDirectory, ...
    [stem '_validation.csv']));

figureHandle = figure('Color','w','Position',[50 50 1430 860]);
layout = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
conditionTitles = {'Mobile: velocity continuity', ...
    'Mobile: no penetration','Mobile: tangential traction jump', ...
    'Cap: relative tangential velocity'};
for conditionIndex = 1:4
    ax = nexttile(layout,conditionIndex);
    semilogy(ax,actualN,max(residualSpeed(:,conditionIndex),eps), ...
        'bo-','MarkerFaceColor','b','LineWidth',1.1);
    hold(ax,'on');
    semilogy(ax,actualN,max(residualForce(:,conditionIndex),eps), ...
        'rs-','MarkerFaceColor','r','LineWidth',1.1);
    title(ax,conditionTitles{conditionIndex});
    ylabel(ax,'Maximum normalized residual');
    xlabel(ax,'N per proxy'); grid(ax,'on');
    if conditionIndex==1
        legend(ax,'Prescribed speed','Prescribed force', ...
            'Location','best');
    end
end
ax = nexttile(layout,5);
semilogy(ax,actualN,max(dragError,eps),'bo-', ...
    'MarkerFaceColor','b','LineWidth',1.1);
title(ax,'Drag error, equation (50)');
ylabel(ax,'|drag/drag_{exact}-1|');
xlabel(ax,'N per proxy'); grid(ax,'on');
ax = nexttile(layout,6);
semilogy(ax,actualN,max(speedError,eps),'rs-', ...
    'MarkerFaceColor','r','LineWidth',1.1);
title(ax,'Speed error, equation (50)');
ylabel(ax,'|U_z/U_{z,exact}-1|');
xlabel(ax,'N per proxy'); grid(ax,'on');
sgtitle(layout,sprintf(['Fixed cap angle %.1f deg; extra collocation %s; ', ...
    'validation minimum edge distance %.3g a'], ...
    rad2deg(options.phi),modeName, ...
    validation.minimumEdgeDistance/options.radius));
exportgraphics(figureHandle,fullfile(outputDirectory,[stem '.png']), ...
    'Resolution',180);
