%% OCT23_STAGNANT_CAP_RING_LEVELS_SWEEP: edge refinement at fixed N and phi
% Investigates whether adding cap-edge collocation rings improves boundary
% accuracy, drag, and gravitational settling speed without changing either
% uniform proxy source set or the base spherical-design collocation nodes.
% Vary only ringLevels: each positive level adds one cap-side and one
% mobile-side ring, increasing M and approaching the singular edge.
% Use one independent validation geometry for both solve modes at every level.
% The eight-panel figure shows four boundary maxima, equation (50) drag
% and speed errors, total M, and closest collocation distance to the edge.
% Edit sourceRequest, phi, and ringLevels below; include level 0 for a
% base-only reference. Figure and CSV files go to experiments/results/.

sourceRequest = 200; % Minimum N on EACH uniform proxy surface.
phi = pi/2;
ringLevels = 1:6;
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
    'sourceRequest',sourceRequest,'collocationRatio',1.5, ...
    'svdTolerance',1e-12,'useExtraCollocation',true, ...
    'gravity',[0 0 -1],'densityDifference',1));
assert(options.phi>0 && options.phi<pi);
assert(isvector(ringLevels) && ~isempty(ringLevels) && ...
    all(ringLevels>=0) && all(ringLevels==round(ringLevels)) && ...
    all(diff(ringLevels)>0));
assert(prescribedSpeed~=0 && isfinite(prescribedSpeed));
ringLevels = ringLevels(:);
Fg = (4*pi*options.radius^3/3)* ...
    options.densityDifference*options.gravity;
exactDrag = stagnant_cap_exact_drag(options.phi, ...
    options.muIn,options.muOut,options.radius);
exactSpeed = Fg(3)/exactDrag;
assert(exactSpeed~=0);

caseCount = numel(ringLevels);
caseOptions = cell(caseCount,1);
geometries = cell(caseCount,1);
for k = 1:caseCount
    o = options;
    o.ringLevels = ringLevels(k);
    caseOptions{k} = stagnant_cap_options(o);
    geometries{k} = stagnant_cap_geometry(caseOptions{k});
    if k>1
        assert(geometries{k}.N==geometries{1}.N && ...
            isequal(geometries{k}.innerSources, ...
                geometries{1}.innerSources) && ...
            isequal(geometries{k}.outerSources, ...
                geometries{1}.outerSources), ...
            'Changing ringLevels must leave both uniform proxy sets fixed.');
    end
end
N = geometries{1}.N;
maxM = max(cellfun(@(g) g.collocationCount,geometries));
% One independent validation geometry serves every level; its finest ring is
% half the finest nonzero collocation offset in this sweep.
finestOffset = options.edgeBandAngle* ...
    2^(-max(1,max(ringLevels)-1));
validation = stagnant_cap_validation_grid(options,maxM,finestOffset);
assert(any(validation.cap) && any(validation.mobile));
fprintf(['Fixed N=%d, phi=%.5g rad. Independent validation: %d points, ', ...
    'minimum edge distance %.6g a.\n'], ...
    N,options.phi,validation.count, ...
    validation.minimumEdgeDistance/options.radius);
fprintf(['Velocity residuals use U_*; traction residuals use ', ...
    'mu_out U_*/a. Maxima include near-edge points.\n']);

collocationM = zeros(caseCount,1);
baseCollocationCount = zeros(caseCount,1);
edgeRingAzimuthCount = zeros(caseCount,1);
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
% Four boundary maxima: mobile continuity, mobile normal velocity,
% mobile tangential traction jump, cap exterior relative tangent velocity.
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
    collocationM(k) = geometry.collocationCount;
    baseCollocationCount(k) = geometry.baseCount;
    edgeRingAzimuthCount(k) = geometry.ringAzimuthCount;
    edgeRingNodeCount(k) = collocationM(k)-baseCollocationCount(k);
    assert(edgeRingNodeCount(k)== ...
        2*ringLevels(k)*edgeRingAzimuthCount(k));
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
    fprintf('\nringLevels=%d: M=%d = %d base + %d ring nodes; nearest %.4g a\n', ...
        ringLevels(k),collocationM(k),baseCollocationCount(k), ...
        edgeRingNodeCount(k),collocationMinimumEdgeDistance(k)/o.radius);
    fprintf('  mobile max [continuity, normal, traction], cap max [tangent]\n');
    fprintf('  speed solve: %.3g  %.3g  %.3g  %.3g\n',residualSpeed(k,:));
    fprintf('  force solve: %.3g  %.3g  %.3g  %.3g\n',residualForce(k,:));
    fprintf('  drag error %.3g, speed error %.3g\n', ...
        dragError(k),speedError(k));
    speedRows = checkSpeed.table;
    forceRows = checkForce.table;
    speedRows.ringLevels = repmat(ringLevels(k),height(speedRows),1);
    speedRows.mode = repmat("speed",height(speedRows),1);
    forceRows.ringLevels = repmat(ringLevels(k),height(forceRows),1);
    forceRows.mode = repmat("force",height(forceRows),1);
    allValidation = [allValidation;speedRows;forceRows]; %#ok<AGROW>
end
if all(ringLevels>=1)
    assert(all(diff(collocationM)>0), ...
        'With fixed N and phi, each extra ring level must increase M.');
end

requestedN = repmat(sourceRequest,caseCount,1);
actualN = repmat(N,caseCount,1);
validationCount = repmat(validation.count,caseCount,1);
validationMinimumEdgeDistance = repmat( ...
    validation.minimumEdgeDistance,caseCount,1);
referenceDrag = repmat(exactDrag,caseCount,1);
referenceSpeed = repmat(exactSpeed,caseCount,1);
summary = table(ringLevels,requestedN,actualN,collocationM, ...
    baseCollocationCount,edgeRingAzimuthCount,edgeRingNodeCount, ...
    collocationMinimumEdgeDistance,validationCount, ...
    validationMinimumEdgeDistance,svdRows,svdCols,svdInputMB, ...
    svdKnownDenseMB,rankSpeed,rankForce,velocityScaleSpeed, ...
    velocityScaleForce,tractionScaleSpeed,tractionScaleForce, ...
    drag,referenceDrag,dragError,speed,referenceSpeed,speedError);
summary.mobileContinuitySpeed = residualSpeed(:,1);
summary.mobileContinuityForce = residualForce(:,1);
summary.mobileNormalSpeed = residualSpeed(:,2);
summary.mobileNormalForce = residualForce(:,2);
summary.mobileTractionSpeed = residualSpeed(:,3);
summary.mobileTractionForce = residualForce(:,3);
summary.capTangentialSpeed = residualSpeed(:,4);
summary.capTangentialForce = residualForce(:,4);
disp(summary(:,{'ringLevels','collocationM','edgeRingNodeCount', ...
    'mobileContinuitySpeed','mobileNormalSpeed', ...
    'mobileTractionSpeed','capTangentialSpeed', ...
    'mobileContinuityForce','mobileNormalForce', ...
    'mobileTractionForce','capTangentialForce', ...
    'dragError','speedError'}));

stem = sprintf('oct23_stagnant_cap_ring_levels_N%04d_phi%05.1f', ...
    N,rad2deg(options.phi));
writetable(summary,fullfile(outputDirectory,[stem '.csv']));
writetable(allValidation,fullfile(outputDirectory, ...
    [stem '_validation.csv']));

figureHandle = figure('Color','w','Position',[50 50 1570 900]);
layout = tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
conditionTitles = {'Mobile: velocity continuity', ...
    'Mobile: no penetration','Mobile: tangential traction jump', ...
    'Cap: relative tangential velocity'};
for conditionIndex = 1:4
    ax = nexttile(layout,conditionIndex);
    semilogy(ax,ringLevels,max(residualSpeed(:,conditionIndex),eps), ...
        'bo-','MarkerFaceColor','b','LineWidth',1.1);
    hold(ax,'on');
    semilogy(ax,ringLevels,max(residualForce(:,conditionIndex),eps), ...
        'rs-','MarkerFaceColor','r','LineWidth',1.1);
    title(ax,conditionTitles{conditionIndex});
    ylabel(ax,'Maximum normalized residual');
    xlabel(ax,'Ring levels'); grid(ax,'on');
    if conditionIndex==1
        legend(ax,'Prescribed speed','Prescribed force', ...
            'Location','best');
    end
end
ax = nexttile(layout,5);
semilogy(ax,ringLevels,max(dragError,eps),'bo-', ...
    'MarkerFaceColor','b','LineWidth',1.1);
title(ax,'Drag error, equation (50)');
ylabel(ax,'|drag/drag_{exact}-1|');
xlabel(ax,'Ring levels'); grid(ax,'on');
ax = nexttile(layout,6);
semilogy(ax,ringLevels,max(speedError,eps),'rs-', ...
    'MarkerFaceColor','r','LineWidth',1.1);
title(ax,'Speed error, equation (50)');
ylabel(ax,'|U_z/U_{z,exact}-1|');
xlabel(ax,'Ring levels'); grid(ax,'on');
ax = nexttile(layout,7);
plot(ax,ringLevels,collocationM,'ko-','MarkerFaceColor','k');
title(ax,'Total collocation count'); ylabel(ax,'M');
xlabel(ax,'Ring levels'); grid(ax,'on');
ax = nexttile(layout,8);
semilogy(ax,ringLevels,collocationMinimumEdgeDistance/options.radius, ...
    'ko-','MarkerFaceColor','k');
title(ax,'Closest collocation point to edge');
ylabel(ax,'min |\theta-\phi| (rad)');
xlabel(ax,'Ring levels'); grid(ax,'on');
sgtitle(layout,sprintf(['Fixed N=%d, phi=%.1f deg, proxy radii %.2ga/%.2ga; ', ...
    'validation minimum edge distance %.3g a'], ...
    N,rad2deg(options.phi),options.innerProxyRatio, ...
    options.outerProxyRatio, ...
    validation.minimumEdgeDistance/options.radius));
exportgraphics(figureHandle,fullfile(outputDirectory,[stem '.png']), ...
    'Resolution',180);
