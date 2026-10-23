%DEMO_MOBILITY_WITH_IMAGES Two-sphere resistance-to-mobility image check.
% Run setpath first, then this demo. Geometry and source/cap controls match
% demo_resistance_with_images.m. The independent check uses more points.
% First compute forces from known rigid motions, then recover those motions
% with mobility. A second discretisation checks the matched-grid result.

clear;
close all;
%% Configure geometry and image discretisation
useSVD = true; % Dense SVD for both resistance and mobility when true.
% Set either flag false to use identity rows in that solver.
useResistanceRowScaling = true;
useMobilityRowScaling = false;
q = [0,0,0;2.001,0,0];
Uref = [1;0;0;0;0;0;-1;0;0;0;0;0];
N = 686;
N_image = []; % If empty: adaptive from image paper Eq. (39).
% imageFamilies tokens: S Stokeslet, R rotlet, D potential dipole, T legacy
% stresslet, T_t line-parallel, T_n first transverse, T_n2 second transverse
% stresslet, Q_s trace-free axial Stokeslet quadrupole. T excludes T_*.
imageFamilies = "SRDQ_s";
capMode = "exponential"; % Switch to "fibonacci" for the paper's caps.
capPointCounts = [];
if capMode == "exponential"
    capPointCounts = 1078;
    rowWeightMode = "local";
else
    rowWeightMode = "paper";
end
opt = struct('N',N,'N_image',N_image, ...
    'imageFamilies',imageFamilies, ... % See token legend above.
    'capMode',capMode,'capPointCounts',capPointCounts, ...
    'rowWeightMode',rowWeightMode, ...
    'capLayers',[18,18],'capExponent',[3,3],'capMinAzimuth',8, ...
    'validationCount',1000,'validationGapCount',401);
% Both solvers and the refined mobility check use this one family selection.
resOpt = opt;
resOpt.useRowScaling = useResistanceRowScaling;
mobOpt = opt;
mobOpt.useRowScaling = useMobilityRowScaling;
if useSVD
    resOpt.solverMode = "denseSVD";
    mobOpt.solverMode = "denseSVD";
end
fprintf('Row scaling: resistance %d, mobility %d\n', ...
    useResistanceRowScaling,useMobilityRowScaling);

%% Compute the resistance solution for the prescribed rigid motions
% F uses the source-strength (force-on-fluid) sign required by mobility.
[F,~,~,~,~] = solve_resistance_with_images(q,Uref,resOpt);

%% Solve mobility and compare recovered motions and force/torque moments
[U,iters,~,~,details] = solve_mobility_with_images(q,F,mobOpt);
translationRows = reshape(6*(0:1)+(1:3).',[],1);
rotationRows = reshape(6*(0:1)+(4:6).',[],1);
if useSVD
    fprintf('Two-sphere mobility: dense SVD rank %d, residual %.3e\n', ...
        details.retainedRanks,details.linearSystemResidual);
else
    fprintf('Two-sphere mobility: GMRES %d, true residual %.3e\n', ...
        iters,details.trueGmresResidual);
end
fprintf('Translation error %.3e; rotation error %.3e\n', ...
    norm(U(translationRows)-Uref(translationRows),inf), ...
    norm(U(rotationRows)-Uref(rotationRows),inf));
fprintf('Maximum force/torque mismatch %.3e\n', ...
    norm(details.realisedForceTorque-F,inf));

%% Check boundary velocity on independent full-sphere and gap-cap points
% Use a different spherical design and a denser Fibonacci gap cap than the fit.
checkOpt = details.options;
fitCounts = arrayfun(@(body) size(body.collocation,1),details.grid.body);
checkOpt.validationCount = max(3000,max(fitCounts)+1000);
checkOpt.validationGapCount = 801;
[err,points,~,relative,gapMask] = image_surface_validation( ...
    q,details.grid,details.lambda,U,checkOpt);
fprintf('Independent surface residual %.3e (%d gap points)\n', ...
    err,sum(cellfun(@nnz,gapMask)));

%% Repeat mobility with a different, refined discretisation
% A different, refined source/collocation discretisation tests the
% resistance-to-mobility answer beyond a matched-grid round trip.
alternative = mobOpt;
alternative.N = 969;
if capMode == "exponential"
    alternative.capMode = "fibonacci";
    alternative.rowWeightMode = "paper";
    alternative.capPointCounts = [];
else
    alternative.capMode = "exponential";
    alternative.rowWeightMode = "local";
    alternative.capPointCounts = 1078;
end
alternative.capLayers = [18,18];
alternative.capExponent = [3,3];
[Ualt,~,~,~,altDetails] = solve_mobility_with_images(q,F,alternative);
[altErr,~,~,~,~] = image_surface_validation( ...
    q,altDetails.grid,altDetails.lambda,Ualt,checkOpt);
fprintf('Alternative grid: translation difference %.3e, rotation difference %.3e\n', ...
    norm(Ualt(translationRows)-U(translationRows),inf), ...
    norm(Ualt(rotationRows)-U(rotationRows),inf));
fprintf('Alternative independent surface residual %.3e\n',altErr);

%% Visualise the independent surface residual, including the gap regions
fig = figure('Name','Two-sphere image mobility: independent residual');
ax = axes(fig);
hold(ax,'on');
for k = 1:2
    x = points{k};
    scatter3(ax,x(:,1),x(:,2),x(:,3),8, ...
        log10(max(relative{k},eps)),'filled');
end
axis(ax,'equal'); grid(ax,'on'); view(ax,3);
xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
colormap(ax,turbo);
cb = colorbar(ax);
cb.Label.String = 'log_{10} relative surface residual';
title(ax,sprintf('P=2, gap=%.2g, independent max residual %.2e', ...
    norm(q(2,:)-q(1,:))-2,err));
