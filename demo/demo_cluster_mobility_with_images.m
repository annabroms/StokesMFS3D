%DEMO_CLUSTER_MOBILITY_WITH_IMAGES Five close spheres in image mobility.
% Geometry, seed, motion, image selection and caps match
% demo_cluster_resistance_with_images.m (which executes P=5).
% Direct summation and body-specific SVDs make this a substantial run.
% Resistance supplies known forces; mobility recovers the original motions.

%% Configure paths, cluster geometry, and image discretisation
clear;
close all;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(root));
P = 5;
gap = 1e-2;
% Set either flag false to use identity rows in that solver.
useResistanceRowScaling = true;
useMobilityRowScaling = false;
capMode = "fibonacci"; % Switch to "exponential" for the graded cap.
rng(20261016);
[q,~] = grow_cluster(P,gap,3);
Uref = randn(6*P,1);
if capMode == "exponential"
    rowWeightMode = "local";
else
    rowWeightMode = "paper";
end
N = 686;
N_image = [];
% imageFamilies tokens: S Stokeslet, R rotlet, D potential dipole, T legacy
% stresslet, T_t line-parallel, T_n first transverse, T_n2 second transverse
% stresslet, Q_s trace-free axial Stokeslet quadrupole. T excludes T_*.
opt = struct('N',N,'proxyRadius',0.63, ...
    'N_image',N_image,'imageFamilies',"SRDT", ...
    'capMode',capMode,'capAngles',[pi/5,pi/60], ...
    'capUpsampling',[6,6],'capLayers',[18,18], ...
    'capExponent',[3,3],'rowWeightMode',rowWeightMode, ...
    'validationCount',180,'validationGapCount',301, ...
    'gmresTolerance',1e-6,'maxIterations',250);
resOpt = opt;
resOpt.useRowScaling = useResistanceRowScaling;
mobOpt = opt;
mobOpt.useRowScaling = useMobilityRowScaling;
fprintf('Row scaling: resistance %d, mobility %d\n', ...
    useResistanceRowScaling,useMobilityRowScaling);

%% Compute the resistance solution for the prescribed rigid motions
% F uses the source-strength (force-on-fluid) sign required by mobility.
[F,~,~,~,~] = solve_resistance_with_images(q,Uref,resOpt);

%% Solve mobility and compare recovered motions and force/torque moments
[U,iters,~,~,details] = solve_mobility_with_images(q,F,mobOpt);
translationRows = reshape(6*(0:P-1)+(1:3).',[],1);
rotationRows = reshape(6*(0:P-1)+(4:6).',[],1);
fprintf('P=%d, gap=%.0e: mobility GMRES %d, residual %.3e\n', ...
    P,gap,iters,details.trueGmresResidual);
fprintf('Translation error %.3e; rotation error %.3e\n', ...
    norm(U(translationRows)-Uref(translationRows),inf), ...
    norm(U(rotationRows)-Uref(rotationRows),inf));
fprintf('Maximum force/torque mismatch %.3e\n', ...
    norm(details.realisedForceTorque-F,inf));

%% Check boundary velocity on independent full-sphere and gap-cap points
% Increase the full-sphere and gap samples beyond the fitted collocation.
checkOpt = details.options;
fitCounts = arrayfun(@(body) size(body.collocation,1),details.grid.body);
checkOpt.validationCount = max(3000,max(fitCounts)+1000);
checkOpt.validationGapCount = 801;
[err,points,~,relative,gapMask] = image_surface_validation( ...
    q,details.grid,details.lambda,U,checkOpt);
fprintf('Independent surface residual %.3e (%d gap points)\n', ...
    err,sum(cellfun(@nnz,gapMask)));

%% Repeat mobility with a different cap discretisation
% Keep the same physical geometry and forces, but change the cap grid.
alternative = mobOpt;
if capMode == "exponential"
    alternative.capMode = "fibonacci";
    alternative.rowWeightMode = "paper";
else
    alternative.capMode = "exponential";
    alternative.rowWeightMode = "local";
end
alternative.capLayers = [18,18];
alternative.capExponent = [3,3];
[Ualt,~,~,~,altDetails] = solve_mobility_with_images(q,F,alternative);
[altErr,~,~,~,~] = image_surface_validation( ...
    q,altDetails.grid,altDetails.lambda,Ualt,checkOpt);
fprintf('Alternative cap grid: translation difference %.3e, rotation difference %.3e\n', ...
    norm(Ualt(translationRows)-U(translationRows),inf), ...
    norm(Ualt(rotationRows)-U(rotationRows),inf));
fprintf('Alternative independent surface residual %.3e\n',altErr);

%% Visualise the independent surface residual, including the gap regions
fig = figure('Name','Five-sphere image mobility: independent residual');
ax = axes(fig);
hold(ax,'on');
for k = 1:P
    x = points{k};
    scatter3(ax,x(:,1),x(:,2),x(:,3),8, ...
        log10(max(relative{k},eps)),'filled');
end
axis(ax,'equal'); grid(ax,'on'); view(ax,3);
xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
colormap(ax,turbo);
cb = colorbar(ax);
cb.Label.String = 'log_{10} relative surface residual';
title(ax,sprintf('P=%d, gap=%.0e, independent max residual %.2e', ...
    P,gap,err));
