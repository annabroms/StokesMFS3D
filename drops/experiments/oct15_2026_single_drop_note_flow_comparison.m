% OCT15_2026_SINGLE_DROP_NOTE_FLOW_COMPARISON
% Compare one spherical-drop MFS flow with Clusters_bisphericals.pdf, p. 2,
% equations (1a)-(1d) and (2). Run this script from any working directory.
% The note uses radius 1, exterior viscosity 1, interior viscosity lambda,
% active stress e_z e_z, and A = 1/[5(lambda+1)]. Its velocity is in the
% laboratory frame and the symmetric drop has zero translational velocity.
% Here the outward normal points from the drop into the exterior. The note's
% tangential jump is (sigma_out-sigma_in)*n = cos(theta)*(e_z-cos(theta)*n),
% so the existing drop solver uses S1=1. Nodal vectors are ordered
% [x1;y1;z1;x2;y2;z2;...].

%% One-drop discretization and independent comparison points
interiorViscosity = 10;
N = 400;
M = 1.2*N;
validationCountRequested = 1000;
innerProxyRadius = 0.4;
outerProxyRadius = 2;
svdTolerance = 1e-12;  % retain sigma > svdTolerance*sigma_max
shellRadii = [0.25 0.5 0.75 0.95 1.05 1.25 1.5 2 3];
planeGridCount = 101;
planeHalfWidth = 2.5;
evaluationBatchSize = 256;
quiverStride = 8;

experimentDirectory = fileparts(mfilename('fullpath'));
dropsDirectory = fileparts(experimentDirectory);
repositoryDirectory = fileparts(dropsDirectory);
addpath(dropsDirectory,fullfile(repositoryDirectory,'common','dense'), ...
    fullfile(repositoryDirectory,'geometry','spheres'));

assert(interiorViscosity > 0 && isfinite(interiorViscosity));
assert(0 < innerProxyRadius && innerProxyRadius < 1 && ...
    1 < outerProxyRadius);
assert(all(shellRadii > 0 & shellRadii ~= 1));
assert(svdTolerance > 0 && svdTolerance < 1);
assert(isscalar(N) && isfinite(N) && N > 0);
assert(isscalar(M) && isfinite(M) && M > 0);
assert(isscalar(validationCountRequested) && ...
    isfinite(validationCountRequested) && validationCountRequested > 0);

% get_sphdesign rounds downward. Select the nearest catalogued count first
% so N, M, and the validation request can be changed independently.
catalogueFile = fullfile(repositoryDirectory,'geometry','spheres', ...
    'spherical_designs','filelist.txt');
catalogueID = fopen(catalogueFile,'r');
assert(catalogueID >= 0,'Could not read the spherical-design catalogue.');
catalogueData = textscan(catalogueID,'%*[^.].%d');
fclose(catalogueID);
availableCounts = catalogueData{1};
assert(~isempty(availableCounts));
[~,sourceIndex] = min(abs(availableCounts-N));
[~,collocationIndex] = min(abs(availableCounts-M));
[~,validationIndex] = min(abs(availableCounts-validationCountRequested));
[sourceDirections,~] = get_sphdesign(availableCounts(sourceIndex));
[surfacePoints,~] = get_sphdesign(availableCounts(collocationIndex));
[validationDirections,~] = get_sphdesign(availableCounts(validationIndex));
sourceCount = size(sourceDirections,1);
collocationCount = size(surfacePoints,1);
validationCount = size(validationDirections,1);
fprintf('Requested/actual spherical-design counts: N %.6g/%d, M %.6g/%d, validation %.6g/%d.\n', ...
    N,sourceCount,M,collocationCount,validationCountRequested,validationCount);

% Rotate the validation directions so none are collocation points.
azimuth = 0.37;
elevation = 0.23;
rotationZ = [cos(azimuth) -sin(azimuth) 0; ...
             sin(azimuth) cos(azimuth) 0; 0 0 1];
rotationY = [cos(elevation) 0 sin(elevation); ...
             0 1 0; -sin(elevation) 0 cos(elevation)];
validationDirections = (rotationZ*rotationY*validationDirections.').';
innerSources = innerProxyRadius*sourceDirections;
outerSources = outerProxyRadius*sourceDirections;

%% Solve the same force-free, fixed-sphere interface problem as the note
% A maps [alpha(3N);gamma(3N)] to [continuity(3M);normal(M);
% tangent traction(2M)]. lambda=alpha-mean(alpha) is force free and
% U=-completionScale*sum(alpha) is the translation completion.
completionScale = 1/(8*pi);
[systemMatrix,tangentOne,tangentTwo] = drop_recompleted_matrix( ...
    innerSources,outerSources,surfacePoints,surfacePoints, ...
    interiorViscosity,completionScale);
cosTheta = surfacePoints(:,3);
tractionData = [cosTheta.*tangentOne(:,3), ...
                cosTheta.*tangentTwo(:,3)];
rhs = [zeros(4*collocationCount,1);reshape(tractionData.',[],1)];

% Radius, forcing, and exterior viscosity are one, so the velocity and
% traction row scales are both one. Keep the physical matrix unchanged.
[leftVectors,singularMatrix,rightVectors] = svd(systemMatrix,'econ');
singularValues = diag(singularMatrix);
retained = singularValues > svdTolerance*singularValues(1);
assert(any(retained));
coefficients = rightVectors(:,retained)* ...
    ((leftVectors(:,retained).'*rhs)./singularValues(retained));
collocationResidual = norm(systemMatrix*coefficients-rhs)/norm(rhs);
alphaNodes = reshape(coefficients(1:3*sourceCount),3,sourceCount);
lambdaNodes = alphaNodes-mean(alphaNodes,2);
lambda = lambdaNodes(:);
gamma = coefficients(3*sourceCount+1:end);
translation = -completionScale*sum(alphaNodes,2).';
noteAmplitude = 1/(5*(interiorViscosity+1));

%% Compare full vector velocities away from the interface
shellMaximumError = zeros(numel(shellRadii),1);
shellRelativeError = zeros(numel(shellRadii),1);
shellReferenceSpeed = zeros(numel(shellRadii),1);
shellPointRelativeError = cell(numel(shellRadii),1);
for shellIndex = 1:numel(shellRadii)
    radius = shellRadii(shellIndex);
    points = radius*validationDirections;
    numericalVelocity = evaluate_drop_field( ...
        points,innerSources,outerSources,lambda,gamma, ...
        interiorViscosity,evaluationBatchSize);
    referenceVelocity = note_drop_velocity(points,noteAmplitude);
    pointError = vecnorm(numericalVelocity-referenceVelocity,2,2);
    shellMaximumError(shellIndex) = max(pointError);
    shellReferenceSpeed(shellIndex) = max(vecnorm(referenceVelocity,2,2));
    shellRelativeError(shellIndex) = ...
        shellMaximumError(shellIndex)/shellReferenceSpeed(shellIndex);
    shellPointRelativeError{shellIndex} = ...
        pointError/shellReferenceSpeed(shellIndex);
end
region = repmat("interior",numel(shellRadii),1);
region(shellRadii > 1) = "exterior";
shellComparison = table(region,shellRadii(:),shellReferenceSpeed, ...
    shellMaximumError,shellRelativeError, ...
    'VariableNames',{'Region','Radius','ReferenceMaxSpeed', ...
    'MaxVectorError','MaxErrorOverReferenceMaxSpeed'});
disp(shellComparison);

% Each panel shows an actual comparison shell. Point color is the full
% vector error divided by that shell's maximum reference speed.
figure('Name','Single drop: comparison shells');
tiledlayout(ceil(numel(shellRadii)/3),3);
for shellIndex = 1:numel(shellRadii)
    nexttile;
    radius = shellRadii(shellIndex);
    points = radius*validationDirections;
    scatter3(points(:,1),points(:,2),points(:,3),12, ...
        shellPointRelativeError{shellIndex},'filled');
    axis equal tight;
    view(35,25);
    grid on;
    xlabel('x'); ylabel('y'); zlabel('z');
    title(sprintf('%s r=%.3g; max %.2e', ...
        region(shellIndex),radius,shellRelativeError(shellIndex)));
    colorbar;
end
sgtitle('Shell samples: |u_{MFS}-u_{note}| / max_{shell}|u_{note}|');

% Check the surface data independently of the flow-field comparison.
surfaceError = drop_validate_surface_conditions( ...
    {lambda},{gamma},translation,{innerSources},{outerSources}, ...
    {validationDirections},{validationDirections},[0 0 1],1, ...
    interiorViscosity,evaluationBatchSize);
noteSurfaceVelocity = note_drop_velocity( ...
    validationDirections,noteAmplitude);
exteriorSurfaceVelocity = reshape(drop_velocity_block( ...
    innerSources,validationDirections,1)*lambda,3,[]).';
interiorSurfaceVelocity = reshape(drop_velocity_block( ...
    outerSources,validationDirections,interiorViscosity)*gamma,3,[]).';
exteriorSurfaceError = max(vecnorm( ...
    exteriorSurfaceVelocity-noteSurfaceVelocity,2,2));
interiorSurfaceError = max(vecnorm( ...
    interiorSurfaceVelocity-noteSurfaceVelocity,2,2));
fprintf('Note amplitude A = %.8g; retained SVD rank %d/%d.\n', ...
    noteAmplitude,nnz(retained),numel(singularValues));
fprintf('Relative collocation residual = %.3e.\n',collocationResidual);
fprintf('Drop translation U = [%.3e %.3e %.3e] (note: zero).\n',translation);
fprintf('Independent surface maxima [continuity, normal, tangent] = [%.3e %.3e %.3e].\n', ...
    surfaceError);
fprintf('Surface velocity error versus note [outside, inside] = [%.3e %.3e].\n', ...
    exteriorSurfaceError,interiorSurfaceError);
fprintf('Maximum shell relative vector error = %.3e.\n',max(shellRelativeError));

%% Compare numerical and analytic flow on the x-z plane
planeAxis = linspace(-planeHalfWidth,planeHalfWidth,planeGridCount);
[xGrid,zGrid] = meshgrid(planeAxis,planeAxis);
planePoints = [xGrid(:),zeros(numel(xGrid),1),zGrid(:)];
numericalPlaneVelocity = evaluate_drop_field( ...
    planePoints,innerSources,outerSources,lambda,gamma, ...
    interiorViscosity,evaluationBatchSize);
referencePlaneVelocity = note_drop_velocity(planePoints,noteAmplitude);
numericalSpeed = reshape(vecnorm(numericalPlaneVelocity,2,2),size(xGrid));
referenceSpeed = reshape(vecnorm(referencePlaneVelocity,2,2),size(xGrid));
vectorError = reshape(vecnorm( ...
    numericalPlaneVelocity-referencePlaneVelocity,2,2),size(xGrid));
commonSpeedMaximum = max([numericalSpeed(:);referenceSpeed(:)]);

figure('Name','Single drop: MFS versus one-cluster note');
tiledlayout(1,3);
circleAngle = linspace(0,2*pi,300);
sampleIndices = 1:quiverStride:planeGridCount;
for panel = 1:3
    nexttile;
    if panel == 1
        field = numericalSpeed;
        velocity = numericalPlaneVelocity;
        caption = 'MFS velocity';
    elseif panel == 2
        field = referenceSpeed;
        velocity = referencePlaneVelocity;
        caption = 'Analytical velocity'; %Note, equations (1)-(2)';
    else
        field = vectorError;
        caption = '|u_{MFS}-u_{analytical}|';
    end
    imagesc(planeAxis,planeAxis,field);
    set(gca,'YDir','normal');
    axis equal tight;
    hold on;
    if panel < 3
        clim([0 commonSpeedMaximum]);
        horizontalVelocity = reshape(velocity(:,1),size(xGrid));
        verticalVelocity = reshape(velocity(:,3),size(xGrid));
        arrowX = horizontalVelocity(sampleIndices,sampleIndices);
        arrowZ = verticalVelocity(sampleIndices,sampleIndices);
        arrowNorm = hypot(arrowX,arrowZ);
        quiver(xGrid(sampleIndices,sampleIndices), ...
            zGrid(sampleIndices,sampleIndices), ...
            arrowX./max(arrowNorm,eps),arrowZ./max(arrowNorm,eps), ...
            0.55,'k');
    end
    plot(cos(circleAngle),sin(circleAngle),'w-','LineWidth',1.2);
    colorbar;
    xlabel('x'); ylabel('z'); title(caption);
end
sgtitle(sprintf('Single drop, viscosity ratio %.3g; color: speed (left, middle), vector error (right)', ...
    interiorViscosity));

function velocity = evaluate_drop_field( ...
    points,innerSources,outerSources,lambda,gamma,viscosity,batchSize)
% Evaluate each MFS field only in its physical fluid domain, in small batches.
velocity = zeros(size(points));
inside = sum(points.^2,2) <= 1;
for domain = 1:2
    if domain == 1
        indices = find(inside);
        sources = outerSources;
        strengths = gamma;
        fluidViscosity = viscosity;
    else
        indices = find(~inside);
        sources = innerSources;
        strengths = lambda;
        fluidViscosity = 1;
    end
    for first = 1:batchSize:numel(indices)
        selected = indices(first:min(first+batchSize-1,numel(indices)));
        block = drop_velocity_block(sources,points(selected,:),fluidViscosity);
        velocity(selected,:) = reshape(block*strengths,3,[]).';
    end
end
end

function velocity = note_drop_velocity(points,amplitude)
% Evaluate the note's spherical fields as Cartesian vectors, including r=0.
radiusSquared = sum(points.^2,2);
z = points(:,3);
inside = radiusSquared <= 1;
velocity = zeros(size(points));
interiorMultiplier = amplitude/2* ...
    (radiusSquared(inside)+2*z(inside).^2-1);
velocity(inside,:) = interiorMultiplier.*points(inside,:);
velocity(inside,3) = velocity(inside,3)-amplitude/2* ...
    (5*radiusSquared(inside)-3).*z(inside);

outside = ~inside;
radius = sqrt(radiusSquared(outside));
exteriorMultiplier = amplitude/2*( ...
    (5./radius.^7-3./radius.^5).*z(outside).^2 + ...
    1./radius.^3-1./radius.^5);
velocity(outside,:) = exteriorMultiplier.*points(outside,:);
velocity(outside,3) = velocity(outside,3)- ...
    amplitude*z(outside)./radius.^5;
end
