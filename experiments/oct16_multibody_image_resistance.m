%OCT16_MULTIBODY_IMAGE_RESISTANCE Reproduce paper triangle/tetrahedron cases.
% Run this script for selected cases. Set runFullSweep below to true for
% the full tetrahedron sweep. Results remain in the workspace.
% Date: 2026-10-16. Execute full convergence experiments deliberately.
% Reference: A. Broms, A. H. Barnett, and A.-K. Tornberg,
% "Accurate close interactions of Stokes spheres using lubrication-adapted
% image systems," Journal of Computational Physics 523 (2025), 113636.
%
% Example 3: three unit spheres at gap 1e-3 translate toward their triangle
% centre. A coarse N=686, Rp=0.63, adaptive-image solve is compared with a fine
% N=1353, Rp=0.7, 30-image solve on independent exterior-plane targets.
% The paper reports maximum relative velocity difference 7.6e-4.
% Example 4: a regular tetrahedron compares enhanced/coarse-only forces
% and torques with a fine MFS reference. Full mode varies gaps and uses 20
% seeded rigid-velocity vectors per gap. All enhanced caps use the paper's
% Fibonacci geometry and area-class row weights.

clear;
close all;
runFullSweep = false; % Set true for the full tetrahedron sweep.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(root));
gap = 1e-3;
side = 2+gap;
triangle = [0,0,0;side,0,0;side/2,sqrt(3)*side/2,0];
triangleCentre = mean(triangle,1);
triangleVelocity = zeros(18,1);
for body = 1:3
    direction = triangleCentre-triangle(body,:);
    direction = direction/norm(direction);
    triangleVelocity(6*(body-1)+(1:3)) = direction.';
end
baseOpt = struct('N',686,'proxyRadius',0.63, ...
    'N_image',[],'imageFamilies',"SRD",'capMode',"fibonacci", ...
    'capAngles',[pi/5,pi/60],'capUpsampling',[6,6], ...
    'rowWeightMode',"paper",'svdCutoff',5e-12, ...
    'gmresTolerance',1e-6,'validationCount',5000, ...
    'validationGapCount',2001);
[~,coarseIterations,coarseNorm,coarseResidual,coarse] = ...
    solve_resistance_with_images(triangle,triangleVelocity,baseOpt);
fineOpt = baseOpt;
fineOpt.N = 1353;
fineOpt.proxyRadius = 0.7;
fineOpt.N_image = 30;
[~,fineIterations,fineNorm,fineResidual,fine] = ...
    solve_resistance_with_images(triangle,triangleVelocity,fineOpt);
planeResolution = [110,80];
if runFullSweep, planeResolution = [250,180]; end
[X,Y] = meshgrid(linspace(-0.5,side+0.5,planeResolution(1)), ...
    linspace(-0.5,sqrt(3)*side/2+0.5,planeResolution(2)));
plane = [X(:),Y(:),zeros(numel(X),1)];
gapResolution = [61,81];
if runFullSweep, gapResolution = [201,161]; end
nearGap = false(size(plane,1),1);
for firstBody = 1:2
    for secondBody = firstBody+1:3
        edge = triangle(secondBody,:)-triangle(firstBody,:);
        edge = edge/norm(edge);
        tangent = [-edge(2),edge(1),0];
        [along,across] = meshgrid( ...
            linspace(-0.04,0.04,gapResolution(1)), ...
            linspace(-0.05,0.05,gapResolution(2)));
        midpoint = (triangle(firstBody,:)+triangle(secondBody,:))/2;
        patch = midpoint+along(:)*edge+across(:)*tangent;
        plane = [plane;patch]; 
        nearGap = [nearGap;true(size(patch,1),1)];
    end
end
outside = true(size(plane,1),1);
for body = 1:3
    outside = outside & vecnorm(plane-triangle(body,:),2,2) > 1;
end
plane = plane(outside,:);
nearGap = nearGap(outside);
coarseFlow = evaluate_flow(coarse,plane,128);
fineFlow = evaluate_flow(fine,plane,128);
fineMagnitude = vecnorm(fineFlow,2,2);
relativeDifference = vecnorm(coarseFlow-fineFlow,2,2)./ ...
    max(fineMagnitude,eps*max(fineMagnitude));
result.triangle = struct('gap',gap,'coarseResidual',coarseResidual, ...
    'fineResidual',fineResidual,'coarseIterations',coarseIterations, ...
    'fineIterations',fineIterations, ...
    'coarseCoefficientNorm',coarseNorm,'fineCoefficientNorm',fineNorm, ...
    'coarseTrueGmresResidual',coarse.trueGmresResidual, ...
    'fineTrueGmresResidual',fine.trueGmresResidual, ...
    'maxPlaneDifference',max(relativeDifference), ...
    'maxGapDifference',max(relativeDifference(nearGap)), ...
    'paperMaxDifference',7.6e-4,'planePoints',size(plane,1));
fprintf('Triangle: max exterior-plane difference %.3e (paper %.3e).\n', ...
    result.triangle.maxPlaneDifference,result.triangle.paperMaxDifference);
fprintf('  near-gap plane difference %.3e on %d targets.\n', ...
    result.triangle.maxGapDifference,nnz(nearGap));
fprintf(['  coarse/fine BC %.3e/%.3e, GMRES %.3e/%.3e, ', ...
    'iterations %d/%d, |lambda| %.3e/%.3e\n'], ...
    coarseResidual,fineResidual,coarse.trueGmresResidual, ...
    fine.trueGmresResidual,coarseIterations,fineIterations, ...
    coarseNorm,fineNorm);
figure('Name','2026-10-16 triangle exterior flow difference');
layout = tiledlayout(1,2);
ax = nexttile(layout);
scatter(ax,plane(:,1),plane(:,2),6, ...
    log10(max(relativeDifference,eps)),'filled');
axis(ax,'equal'); colorbar(ax); xlabel(ax,'x'); ylabel(ax,'y');
title(ax,'log10 relative coarse/fine exterior velocity difference');
ax = nexttile(layout);
scatter(ax,plane(nearGap,1),plane(nearGap,2),7, ...
    log10(max(relativeDifference(nearGap),eps)),'filled');
axis(ax,'equal'); colorbar(ax);
xlim(ax,[0.93,1.07]); ylim(ax,[-0.06,0.06]);
xlabel(ax,'x'); ylabel(ax,'y'); title(ax,'Close-up of one gap');

tetraGaps = 1e-3;
realizations = 1;
if runFullSweep
    tetraGaps = logspace(-3,-0.3,12);
    realizations = 20;
end
rng(20261016);
result.tetrahedron = repmat(struct('gap',0,'realization',0, ...
    'imageForceError',0,'basicForceError',0, ...
    'imageResidual',0,'basicResidual',0, ...
    'imageIterations',0,'basicIterations',0, ...
    'imageTrueGmresResidual',0,'basicTrueGmresResidual',0, ...
    'imageCoefficientNorm',0,'basicCoefficientNorm',0),0,1);
for delta = tetraGaps
    sideLength = 2+delta;
    tetra = [0,0,0;sideLength,0,0; ...
        sideLength/2,sqrt(3)*sideLength/2,0; ...
        sideLength/2,sideLength/(2*sqrt(3)),sqrt(2/3)*sideLength];
    for realization = 1:realizations
        velocity = randn(24,1);
        enhancedOpt = baseOpt;
        enhancedOpt.N_image = [];
        [imageForce,imageIt,imageNorm,imageRes,imageDetail] = ...
            solve_resistance_with_images(tetra,velocity,enhancedOpt);
        basicOpt = baseOpt;
        basicOpt.N_image = 0;
        [basicForce,basicIt,basicNorm,basicRes,basicDetail] = ...
            solve_resistance_with_images(tetra,velocity,basicOpt);
        referenceOpt = fineOpt;
        [referenceForce,~,~,~] = ...
            solve_resistance_with_images(tetra,velocity,referenceOpt);
        forceScale = max(norm(referenceForce,inf),eps);
        row = struct('gap',delta,'realization',realization, ...
            'imageForceError',norm(imageForce-referenceForce,inf)/forceScale, ...
            'basicForceError',norm(basicForce-referenceForce,inf)/forceScale, ...
            'imageResidual',imageRes,'basicResidual',basicRes, ...
            'imageIterations',imageIt,'basicIterations',basicIt, ...
            'imageTrueGmresResidual',imageDetail.trueGmresResidual, ...
            'basicTrueGmresResidual',basicDetail.trueGmresResidual, ...
            'imageCoefficientNorm',imageNorm, ...
            'basicCoefficientNorm',basicNorm);
        result.tetrahedron(end+1,1) = row;
        fprintf(['Tetra gap %.2e run %d: image/basic force err %.3e/%.3e, ', ...
            'BC %.3e/%.3e, GMRES %.3e/%.3e, iterations %d/%d, ', ...
            '|lambda| %.3e/%.3e\n'],delta,realization, ...
            row.imageForceError,row.basicForceError,imageRes,basicRes, ...
            imageDetail.trueGmresResidual,basicDetail.trueGmresResidual, ...
            imageIt,basicIt,imageNorm,basicNorm);
    end
end
figure('Name','2026-10-16 tetrahedron force and torque errors');
ax = axes; hold(ax,'on'); grid(ax,'on');
gapIndex = 0;
for gapValue = unique([result.tetrahedron.gap])
    gapIndex = gapIndex+1;
    rows = [result.tetrahedron.gap] == gapValue;
    enhancedErrors = [result.tetrahedron(rows).imageForceError];
    basicErrors = [result.tetrahedron(rows).basicForceError];
    imageHandle = loglog(ax,gapValue,median(enhancedErrors), ...
        'bo','MarkerFaceColor','b');
    basicHandle = loglog(ax,gapValue,median(basicErrors), ...
        'rs','MarkerFaceColor','r');
    if gapIndex == 1
        legendHandles = [imageHandle,basicHandle];
    end
    if numel(enhancedErrors) > 1
        loglog(ax,[gapValue,gapValue], ...
            [min(enhancedErrors),max(enhancedErrors)],'b-');
        loglog(ax,[gapValue,gapValue], ...
            [min(basicErrors),max(basicErrors)],'r-');
    end
end
yline(ax,1e-3,'k--','10^{-3} target');
xlabel(ax,'Sphere gap'); ylabel(ax,'Relative force/torque error');
title(ax,'Tetrahedron: medians and ranges against the fine MFS grid');
legend(ax,legendHandles,{'With images','Proxy only'},'Location','best');

function velocity = evaluate_flow(solution,targets,batchSize)
%EVALUATE_FLOW evaluates all stored source families on exterior targets.
velocity = zeros(size(targets));
for body = 1:numel(solution.grid.body)
    for first = 1:batchSize:size(targets,1)
        nodes = first:min(first+batchSize-1,size(targets,1));
        G = stokes_mfs_velocity_matrix( ...
            solution.grid.body(body).basis,targets(nodes,:));
        velocity(nodes,:) = velocity(nodes,:)+ ...
            reshape(G*solution.lambda{body},3,[]).';
    end
end
end
