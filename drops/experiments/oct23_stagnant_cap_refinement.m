%% OCT23_STAGNANT_CAP_REFINEMENT: fixed-physics source-ring comparison.
% Run explicitly from MATLAB. The two cases share physical parameters,
% uniform proxy nodes, SVD tolerance, and an independent validation geometry.
% Edit n_level for cone source rings and ringLevels for added collocation
% rings on each side of the cap edge. The closest ring offset is controlled
% by collocAngularExtent and collocClusterRatio.
% Stagnant-cap refinement has been validated with S (Stokeslets). The base
% proxies are always S. refinementFamilies also accepts R (rotlet),
% D (potential dipole), T (legacy stresslet), T_t/T_n/T_n2 (oriented
% stresslets), and Q_s (axial Stokeslet quadrupole) through the shared image
% kernels. These other choices are experimental for this stagnant-cap solve;
% check independent boundary residuals and the theoretical speed before
% treating their results as accurate. T cannot be combined with T_t/T_n/T_n2.
% R can create source torque, which this nonrotating model reports but does
% not constrain. The selection below uses S and the legacy T stresslet.

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
fprintf('Options: %s\nSource-family parser: %s\n', ...
    which('stagnant_cap_options'),which('image_family_list'));

opt = stagnant_cap_options(struct( ...
    'radius',1,'muOut',1,'muIn',2,'phi',pi/4, ...
    'innerProxyRatio',0.5,'outerProxyRatio',1.8, ...
    'sourceRequest',400,'collocationRatio',1.5, ...
    'svdTolerance',1e-12,'gravity',[0 0 -1], ...
    'densityDifference',1,'useExtraCollocation',false, ...
    'n_level',10,'ringLevels',80, ...
    'innerSourceStartDistance',0.35, ...
    'outerSourceStartDistance',0.45, ...
    'sourceClusterRatio',0.7,'n_az_source',64, ...
    'n_az_colloc',64,'collocAngularExtent',0.2, ...
    'collocClusterRatio',0.8, ...
    'refinementFamilies',["S","T"])); %T_t","T_n","T_n2"]));
% Convert the requested per-region ring count for the shared cone builder.
familyCount = numel(opt.refinementFamilies)+ ...
    double(opt.derivative_direction~="none");
assert(opt.n_level>0 && familyCount>0 && ...
    opt.ringLevels>=2*opt.n_level*familyCount);
opt.colloc_factor = opt.ringLevels/(opt.n_level*familyCount);
fprintf('Enhancing source families: %s\n', ...
    strjoin(opt.refinementFamilies,', '));
Fg = (4*pi*opt.radius^3/3)*opt.densityDifference*opt.gravity;
exactDrag = stagnant_cap_exact_drag( ...
    opt.phi,opt.muIn,opt.muOut,opt.radius);
exactSpeed = Fg(3)/exactDrag;
base = stagnant_cap_geometry(opt);
refined = stagnant_cap_refined_geometry(opt);
assert(isequal(refined.collocationRingCount, ...
    [opt.ringLevels opt.ringLevels]));
assert(isequal(base.innerSources,refined.innerBaseSources) && ...
    isequal(base.outerSources,refined.outerBaseSources) && ...
    isequal(base.points,refined.points(1:refined.baseCount,:)));
% Report the dense solve's memory demand before assembling its matrices.
plannedRows = 6*refined.collocationCount;
plannedColumns = refined.exteriorCoefficientCount + ...
    refined.interiorCoefficientCount;
thinSize = min(plannedRows,plannedColumns);
inputGB = 8*plannedRows*plannedColumns/1e9;
factorGB = 8*(plannedRows*thinSize + ...
    plannedColumns*thinSize + thinSize^2)/1e9;
fprintf(['Refined solve preflight: %d x %d; one matrix %.2f GB; ', ...
    'known dense SVD storage >= %.2f GB, excluding assembly arrays ', ...
    'and workspace.\n'],plannedRows,plannedColumns,inputGB, ...
    2*inputGB+factorGB);
% assert(isequal(refined.addedSourceCount,[320 320]) && ...
%     isequal(refined.addedCollocationCount,[640 640]));
grids = {base,refined};
names = ["uniform";"refined"];

% One shifted, denser validation set is reused without changing its sample
% measure between the two solves. Its closest ring is nearer than either
% refined collocation geometry's closest ring.
validation = stagnant_cap_validation_grid( ...
    opt,refined.collocationCount,refined.ringOffsets(end));
assert(validation.count>=5*refined.collocationCount);
fprintf(['Independent validation (%s base): %d points, ', ...
    'minimum edge distance %.6g a.\n'], ...
    validation.baseMethod,validation.count, ...
    validation.minimumEdgeDistance/opt.radius);
solutions = cell(2,1);
reports = cell(2,1);
validationRows = table();
for k = 1:2
    geometry = grids{k};
    solutions{k} = solve_stagnant_cap( ...
        opt,struct('kind','force','Fg',Fg),geometry);
    reports{k} = stagnant_cap_validate(solutions{k},validation,opt);
    residualFigure = plot_stagnant_cap_residuals( ...
        reports{k},validation,solutions{k},opt,'on');
    exportgraphics(residualFigure,fullfile(outputDirectory, ...
        sprintf('oct23_stagnant_cap_refinement_residuals_%s.png', ...
        names(k))),'Resolution',150);
    rows = reports{k}.table;
    rows.caseName = repmat(names(k),height(rows),1);
    validationRows = [validationRows;rows]; 
end

sourceRingsInside = [0;refined.sourceRingCount(1)];
sourceRingsOutside = [0;refined.sourceRingCount(2)];
ringLevels = [0;opt.ringLevels];
addedSourcesInside = [0;refined.addedSourceCount(1)];
addedSourcesOutside = [0;refined.addedSourceCount(2)];
collocationRingsCap = [0;refined.collocationRingCount(1)];
collocationRingsMobile = [0;refined.collocationRingCount(2)];
addedCollocationCap = [0;refined.addedCollocationCount(1)];
addedCollocationMobile = [0;refined.addedCollocationCount(2)];
sourcePointsInside = zeros(2,1);
sourcePointsOutside = zeros(2,1);
totalSourcePoints = zeros(2,1);
collocationPoints = zeros(2,1);
unknownCount = zeros(2,1);
svdRows = zeros(2,1);
svdColumns = zeros(2,1);
retainedRank = zeros(2,1);
translationZ = zeros(2,1);
speedRelativeError = zeros(2,1);
forceRelativeResidual = zeros(2,1);
torqueRelativeMagnitude = zeros(2,1);
boundaryMaxima = zeros(2,4);
for k = 1:2
    geometry = grids{k};
    solution = solutions{k};
    sourcePointsInside(k) = size(geometry.innerSources,1);
    sourcePointsOutside(k) = size(geometry.outerSources,1);
    totalSourcePoints(k) = ...
        sourcePointsInside(k)+sourcePointsOutside(k);
    collocationPoints(k) = geometry.collocationCount;
    unknownCount(k) = numel(solution.lambda)+numel(solution.gamma);
    svdRows(k) = solution.matrixSize(1);
    svdColumns(k) = solution.matrixSize(2);
    retainedRank(k) = solution.retainedRank;
    translationZ(k) = solution.U(3);
    speedRelativeError(k) = abs(solution.U(3)/exactSpeed-1);
    forceRelativeResidual(k) = solution.forceResidual/norm(Fg);
    torqueRelativeMagnitude(k) = ...
        norm(solution.torque)/(norm(Fg)*opt.radius);
    % The four imposed conditions are evaluated on their full domains.
    boundaryMaxima(k,:) = [ ...
        max(reports{k}.point.continuity)/solution.velocityScale, ...
        max(reports{k}.point.normal)/solution.velocityScale, ...
        max(reports{k}.point.capOuter(validation.cap))/ ...
            solution.velocityScale, ...
        max(reports{k}.point.mobileTraction(validation.mobile))/ ...
            solution.tractionScale];
end
summary = table(names,sourceRingsInside,sourceRingsOutside,ringLevels, ...
    addedSourcesInside,addedSourcesOutside,sourcePointsInside, ...
    sourcePointsOutside,totalSourcePoints,collocationRingsCap, ...
    collocationRingsMobile,addedCollocationCap, ...
    addedCollocationMobile,collocationPoints,unknownCount, ...
    svdRows,svdColumns,retainedRank,translationZ, ...
    speedRelativeError,forceRelativeResidual,torqueRelativeMagnitude);
summary.continuityMaxRelative = boundaryMaxima(:,1);
summary.normalMaxRelative = boundaryMaxima(:,2);
summary.capTangentialMaxRelative = boundaryMaxima(:,3);
summary.mobileTractionMaxRelative = boundaryMaxima(:,4);
summary.allBoundaryMaxRelative = max(boundaryMaxima,[],2);
disp(summary);
fprintf(['Maximum relative boundary errors on the independent surface: ', ...
    'velocity residuals / U_*, traction residual / ', ...
    '(muOut*U_*/a).\n']);
disp(summary(:,{'names','continuityMaxRelative', ...
    'normalMaxRelative','capTangentialMaxRelative', ...
    'mobileTractionMaxRelative','allBoundaryMaxRelative'}));
disp(validationRows(:,{'caseName','condition','region','zone', ...
    'nodeCount','minimumEdgeDistance','maximum','rms'}));
fprintf('Theoretical axial speed %.12g.\n',exactSpeed);
fprintf('Validation minimum edge distance %.6g a.\n', ...
    validation.minimumEdgeDistance/opt.radius);
fprintf('Closest source distance %.6g a; source azimuthal arc %.6g a outside.\n', ...
    refined.sourceDistances.outer(end)/opt.radius, ...
    2*pi*(opt.radius+refined.sourceDistances.outer(end))* ...
    sin(opt.phi)/(opt.n_az_source*opt.radius));
writetable(summary,fullfile(outputDirectory, ...
    'oct23_stagnant_cap_refinement_summary.csv'));
writetable(validationRows,fullfile(outputDirectory, ...
    'oct23_stagnant_cap_refinement_validation.csv'));

% Evaluate both traces at the same independent meridional locations.
near = logspace(-3,log10(opt.collocAngularExtent),90);
theta = sort(unique([linspace(0.02,pi-0.02,180), ...
    opt.phi-near,opt.phi+near]));
theta = theta(abs(theta-opt.phi)>1e-12);
theta = theta(:);
normal = [sin(theta),zeros(numel(theta),1),cos(theta)];
profile.points = opt.centre+opt.radius*normal;
profile.normals = normal;
profile.theta = theta;
profile.cap = theta<opt.phi;
profile.mobile = theta>opt.phi;
profile.edgeDistance = opt.radius*abs(theta-opt.phi);
profile.nearEdge = profile.edgeDistance<= ...
    opt.radius*opt.edgeBandAngle;
profile.minimumEdgeDistance = min(profile.edgeDistance);
profile.count = numel(theta);
velocity = zeros(profile.count,2);
for k = 1:2
    trace = stagnant_cap_validate(solutions{k},profile,opt);
    velocity(:,k) = trace.point.meridionalVelocity;
end
s = opt.radius*(theta(profile.mobile)-opt.phi);
profileFigure = figure('Name','Stagnant-cap refinement velocity', ...
    'Color','w');
tiledlayout(profileFigure,1,3,'TileSpacing','compact');
ax = nexttile; hold(ax,'on');
plot(ax,theta,velocity(:,1),'b-',theta,velocity(:,2),'r-');
xline(ax,opt.phi,'k--');
xlabel(ax,'polar angle \theta'); ylabel(ax,'drop-relative u_\theta');
legend(ax,names,'Location','best'); grid on
ax = nexttile; hold(ax,'on');
plot(ax,s,velocity(profile.mobile,1),'b-', ...
    s,velocity(profile.mobile,2),'r-');
plot(ax,0,0,'ko','MarkerFaceColor','k');
xlabel(ax,'mobile edge distance s'); ylabel(ax,'drop-relative u_\theta');
grid on
ax = nexttile; hold(ax,'on');
semilogx(ax,s,velocity(profile.mobile,1)./sqrt(s),'b-', ...
    s,velocity(profile.mobile,2)./sqrt(s),'r-');
xlabel(ax,'mobile edge distance s'); ylabel(ax,'u_\theta / \surd s');
grid on
exportgraphics(profileFigure,fullfile(outputDirectory, ...
    'oct23_stagnant_cap_refinement_velocity.png'),'Resolution',180);

% Added families share source coordinates, so each needs its own panel.
for k = 1:2
    strengthFigure = plot_stagnant_cap_source_strengths( ...
        solutions{k},names(k),'on');
    exportgraphics(strengthFigure,fullfile(outputDirectory, ...
        sprintf('oct23_stagnant_cap_refinement_strengths_%s.png', ...
        names(k))),'Resolution',180);
end
