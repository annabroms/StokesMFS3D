%% OCT23_STAGNANT_CAP_ANGLE_SWEEP: drag and settling speed versus cap angle
% Investigates how the axial drag and gravitational settling speed vary
% from the clean drop (phi=0) to the fully immobilized sphere (phi=pi).
% Vary phi while holding the physical parameters, proxy radii, and requested
% source count N fixed. The collocation count may change with phi because
% the cap-edge circumference sets the added ring population.
% angleOptions.ringLevels counts added collocation rings per region for
% nondegenerate cap angles; neither endpoint has edge rings.
% At each angle, solve a prescribed-speed resistance problem and a separate
% prescribed-force mobility problem; compare each result with note equation
% (50). The four-panel figure shows both quantities and their relative
% differences; the CSV records counts, SVD sizes, ranks, and errors.
% Edit angleSourceRequest, angleCount, angleOptions, or the collocation
% toggle below. Results are saved under experiments/results/.
clear; close all;
useExtraCollocation = true; % Set false for base spherical-design nodes only.
experimentDirectory = fileparts(mfilename('fullpath'));
repositoryDirectory = fileparts(fileparts(experimentDirectory));
addpath(fullfile(repositoryDirectory,'drops'), ...
    fullfile(repositoryDirectory,'images'), ...
    fullfile(repositoryDirectory,'common','dense'), ...
    fullfile(repositoryDirectory,'geometry','spheres'), ...
    fullfile(repositoryDirectory,'tools'));
outputDirectory = fullfile(experimentDirectory,'results');
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
if ~exist('angleSourceRequest','var')
    angleSourceRequest = 700; 
end
if ~exist('angleCount','var'), angleCount = 13; end
if ~exist('angleOptions','var')
    angleOptions = stagnant_cap_options(struct( ...
        'radius',1,'muOut',1,'muIn',2, ...
        'innerProxyRatio',0.63,'outerProxyRatio',2, ...
        'svdTolerance',1e-12,'gravity',[0 0 -1], ...
        'densityDifference',1,'ringLevels',4, ...
        'collocationRatio',1.5));
end
assert(angleCount>=3 && angleCount==round(angleCount));
angleOptions.useExtraCollocation = useExtraCollocation;
angleOptions = stagnant_cap_options(angleOptions);
gravityForce = (4*pi*angleOptions.radius^3/3)* ...
    angleOptions.densityDifference*angleOptions.gravity;
assert(gravityForce(3)~=0);
angles = linspace(0,pi,angleCount).';
degrees = rad2deg(angles);
actualN = zeros(angleCount,1);
ringLevels = zeros(angleCount,1);
collocationM = zeros(angleCount,1);
baseCollocationCount = zeros(angleCount,1);
edgeRingAzimuthCount = zeros(angleCount,1);
edgeRingNodeCount = zeros(angleCount,1);
svdRows = zeros(angleCount,1);
svdCols = zeros(angleCount,1);
svdInputMB = zeros(angleCount,1);
svdKnownDenseMB = zeros(angleCount,1);
numericalForce = zeros(angleCount,1);
referenceForce = stagnant_cap_exact_drag(angles, ...
    angleOptions.muIn,angleOptions.muOut,angleOptions.radius);
numericalSpeed = zeros(angleCount,1);
referenceSpeed = gravityForce(3)./referenceForce;
rankResistance = zeros(angleCount,1);
rankMobility = zeros(angleCount,1);

for k = 1:angleCount
    o = angleOptions;
    o.phi = angles(k);
    o.sourceRequest = angleSourceRequest;
    o.edgeBandAngle = []; % Recompute the admissible band at this angle.
    o = stagnant_cap_options(o);
    geometry = stagnant_cap_geometry(o);
    resistance = solve_stagnant_cap(o, ...
        struct('kind','speed','Uz',-1),geometry);
    mobility = solve_stagnant_cap(o, ...
        struct('kind','force','Fg',gravityForce),geometry);
    numericalForce(k) = resistance.FH(3);
    numericalSpeed(k) = mobility.U(3);
    actualN(k) = geometry.N;
    ringLevels(k) = numel(geometry.ringOffsets);
    collocationM(k) = geometry.collocationCount;
    baseCollocationCount(k) = geometry.baseCount;
    edgeRingAzimuthCount(k) = geometry.ringAzimuthCount;
    edgeRingNodeCount(k) = collocationM(k)-baseCollocationCount(k);
    svdRows(k) = resistance.matrixSize(1);
    svdCols(k) = resistance.matrixSize(2);
    svdInputMB(k) = resistance.svdInputMB;
    svdKnownDenseMB(k) = resistance.svdKnownDenseMB;
    assert(isequal(resistance.matrixSize,mobility.matrixSize));
    rankResistance(k) = resistance.retainedRank;
    rankMobility(k) = mobility.retainedRank;
    fprintf(['phi=%6.1f deg, N=%d, M=%d (base %d + edge rings %d): ', ...
        'F=%.7g (exact %.7g), U=%.7g (exact %.7g)\n'], ...
        degrees(k),actualN(k),collocationM(k), ...
        baseCollocationCount(k),edgeRingNodeCount(k),numericalForce(k), ...
        referenceForce(k),numericalSpeed(k),referenceSpeed(k));
end
relativeForceDifference = abs(numericalForce-referenceForce)./abs(referenceForce);
relativeSpeedDifference = abs(numericalSpeed-referenceSpeed)./abs(referenceSpeed);
results = table(angles,degrees,actualN,ringLevels,collocationM, ...
    baseCollocationCount,edgeRingAzimuthCount,edgeRingNodeCount, ...
    svdRows,svdCols, ...
    svdInputMB,svdKnownDenseMB,numericalForce, ...
    referenceForce,relativeForceDifference,numericalSpeed,referenceSpeed, ...
    relativeSpeedDifference,rankResistance,rankMobility);
writetable(results,fullfile(outputDirectory, ...
    'oct23_stagnant_cap_angle_sweep.csv'));

fineAngles = linspace(0,pi,401);
fineDegrees = rad2deg(fineAngles);
fineResistance = stagnant_cap_exact_drag(fineAngles, ...
    angleOptions.muIn,angleOptions.muOut,angleOptions.radius);
figureHandle = figure('Color','w','Position',[50 50 1250 820]);
layout = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(layout,1);
plot(ax,fineDegrees,fineResistance,'k-','LineWidth',1.7);
hold(ax,'on');
plot(ax,degrees,numericalForce,'bo','MarkerFaceColor','b');
title(ax,'Drag force for U_z=-1'); ylabel(ax,'-F_H,z/U_z');
legend(ax,'Equation (50)','MFS','Location','best'); grid on;
ax = nexttile(layout,2);
semilogy(ax,degrees,max(relativeForceDifference,eps),'bo-', ...
    'MarkerFaceColor','b');
title(ax,'Drag force relative difference');
ylabel(ax,'|F_{MFS}-F_{exact}|/|F_{exact}|'); grid on;
ax = nexttile(layout,3);
plot(ax,fineDegrees,gravityForce(3)./fineResistance,'k-', ...
    'LineWidth',1.7);
hold(ax,'on');
plot(ax,degrees,numericalSpeed,'ro','MarkerFaceColor','r');
title(ax,sprintf('Speed for F_{g,z}=%.4g',gravityForce(3)));
ylabel(ax,'U_z');
legend(ax,'Equation (50): F_g/drag','MFS','Location','best');
grid on
ax = nexttile(layout,4);
semilogy(ax,degrees,max(relativeSpeedDifference,eps),'ro-', ...
    'MarkerFaceColor','r');
title(ax,'Speed relative difference');
ylabel(ax,'|U_{MFS}-U_{exact}|/|U_{exact}|'); grid on;
for k = 1:4
    ax = nexttile(layout,k);
    xlabel(ax,'Cap angle \phi (degrees)');
    xlim(ax,[0 180]);
end
sgtitle(layout,sprintf('Stagnant cap: mu_{in}/mu_{out}=%.3g, N=%d, proxy radii %.2ga and %.2ga', ...
    angleOptions.muIn/angleOptions.muOut,actualN(1), ...
    angleOptions.innerProxyRatio,angleOptions.outerProxyRatio));
filename = sprintf('oct23_stagnant_cap_angle_sweep_N%u.png',actualN(1));
exportgraphics(figureHandle,fullfile(outputDirectory, ...
    filename),'Resolution',180);
disp(results);
