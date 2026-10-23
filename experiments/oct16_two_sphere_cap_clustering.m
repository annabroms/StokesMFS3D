%OCT16_TWO_SPHERE_CAP_CLUSTERING Isolate collocation-cap accuracy effects.
% Run this script for the matched-cap comparison. Set runFullSweep below to
% true for the full parameter sweep. Results remain in the workspace.
% Date: 2026-10-16. This experiment is intended for deliberate execution.
% Reference: A. Broms, A. H. Barnett, and A.-K. Tornberg,
% "Accurate close interactions of Stokes spheres using lubrication-adapted
% image systems," Journal of Computational Physics 523 (2025), 113636.
%
% Two unit spheres have gap 1e-3 and squeeze along the line of centres.
% One proxy-only solve uses the basic source and collocation designs. All
% Fibonacci and exponential cap solves share the same 20-node S/R/D image
% lines and coarse design; only cap placement and its local row areas vary.
% solverMode selects one-body-preconditioned GMRES or dense weighted SVD.
% A 6000-node independent spherical design plus 3001 gap-cap nodes per body
% checks the boundary condition in every case. Figures show collocation
% locations and log10 relative residuals over the sphere and near the gap.

clear;
close all;
runFullSweep = false; % Set true for the full cap-parameter sweep.
solverMode = "oneBody"; % Choose "oneBody" or "denseSVD".
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(root));
gap = 1e-3;
q = [0,0,0;2+gap,0,0];
U = [1;0;0;0;0;0;-1;0;0;0;0;0];
angleScales = 1;
countScales = 1;
exponents = 3;
if runFullSweep
    angleScales = [0.75,1,1.25];
    countScales = [0.75,1,1.25];
    exponents = [2,3,4];
end
resultTemplate = struct('mode',"",'solverMode',solverMode, ...
    'useImages',false, ...
    'angleScale',0,'countScale',0,'exponent',0, ...
    'maxResidual',0,'gapResidual',0,'iterations',0, ...
    'trueGmresResidual',0,'linearSystemResidual',0, ...
    'coefficientNorm',0,'force',zeros(12,1));
caseCount = 1+numel(angleScales)*numel(countScales)*(1+numel(exponents));
results = repmat(resultTemplate,caseCount,1);
resultIndex = 1;
representativeDetails = cell(1,3);

% The single no-image reference has no cap collocation enhancement.
basicOpt = struct('N',686,'proxyRadius',0.63, ...
    'N_image',0,'useImages',false,'rowWeightMode',"paper", ...
    'solverMode',solverMode, ...
    'validationCount',6000,'validationGapCount',3001);
[F,it,lambdaNorm,err,detail] = ...
    solve_resistance_with_images(q,U,basicOpt);
assert(all(arrayfun(@(body) isempty(body.capPoints),detail.grid.body)));
assert(all(arrayfun(@(body) isempty(body.imagePoints),detail.grid.body)));
near = detail.validationGap{1};
assert(any(near),'The basic case must still sample the gap independently.');
results(resultIndex) = struct('mode',"basic",'solverMode',solverMode, ...
    'useImages',false, ...
    'angleScale',NaN,'countScale',NaN,'exponent',NaN, ...
    'maxResidual',err, ...
    'gapResidual',max(detail.validationRelative{1}(near)), ...
    'iterations',it,'trueGmresResidual',detail.trueGmresResidual, ...
    'linearSystemResidual',detail.linearSystemResidual, ...
    'coefficientNorm',lambdaNorm,'force',F);
representativeDetails{1} = detail;
fprintf(['Proxy only, %s: BC %.3e, relres %.3e, linear residual %.3e, ', ...
    'it %d, |lambda| %.3e\n'],solverMode,err,results(1).gapResidual, ...
    detail.linearSystemResidual,it,lambdaNorm);

% Hold the image representation fixed while comparing cap discretisations.
modes = ["fibonacci","exponential"];
for modeIndex = 1:2
    mode = modes(modeIndex);
    modeExponents = exponents;
    if mode == "fibonacci"
        modeExponents = 3;
    end
    % Retain a plotted case even when the requested exponent omits 3.
    representativeExponent = modeExponents(1);
    if any(modeExponents == 3)
        representativeExponent = 3;
    end
    for angleScale = angleScales
        for countScale = countScales
            for exponent = modeExponents
                opt = basicOpt;
                opt.N_image = 20;
                opt.useImages = true;
                % imageFamilies tokens: S Stokeslet, R rotlet, D potential
                % dipole, T legacy stresslet, T_t line-parallel, T_n first
                % transverse, T_n2 second transverse stresslet, Q_s trace-free
                % axial Stokeslet quadrupole. T excludes T_*.
                opt.imageFamilies = "SRD";
                opt.capMode = mode;
                opt.capAngles = angleScale*[pi/5,pi/60];
                opt.capUpsampling = countScale*[6,6];
                opt.capLayers = [18,18];
                opt.capExponent = [exponent,exponent];
                opt.rowWeightMode = "local";
                [F,it,lambdaNorm,err,detail] = ...
                    solve_resistance_with_images(q,U,opt);
                near = detail.validationGap{1};
                entry = struct('mode',mode,'solverMode',solverMode, ...
                    'useImages',true, ...
                    'angleScale',angleScale,'countScale',countScale, ...
                    'exponent',exponent,'maxResidual',err, ...
                    'gapResidual',max(detail.validationRelative{1}(near)), ...
                    'iterations',it, ...
                    'trueGmresResidual',detail.trueGmresResidual, ...
                    'linearSystemResidual',detail.linearSystemResidual, ...
                    'coefficientNorm',lambdaNorm,'force',F);
                resultIndex = resultIndex+1;
                results(resultIndex) = entry;
                fprintf('%s, %s, angle=%g, count=%g, exp=%g: ', ...
                    solverMode,mode,angleScale,countScale,exponent);
                fprintf(['BC %.3e, gap %.3e, linear residual %.3e, ', ...
                    'it %d, |lambda| %.3e\n'],err,entry.gapResidual, ...
                    detail.linearSystemResidual,it,lambdaNorm);
                if angleScale == 1 && countScale == 1 && ...
                        exponent == representativeExponent
                    representativeDetails{modeIndex+1} = detail;
                end
            end
        end
    end
end

labels = ["Proxy only","Fibonacci images","Exponential images"];
residualFigure = figure('Name','2026-10-16 cap study: surface residual');
residualLayout = tiledlayout(residualFigure,2,3,'TileSpacing','compact');
distributionFigure = figure('Name','2026-10-16 cap study: collocation');
distributionLayout = tiledlayout(distributionFigure,1,3,'TileSpacing','compact');
closeFigure = figure('Name','2026-10-16 cap study: contact collocation');
closeLayout = tiledlayout(closeFigure,1,3,'TileSpacing','compact');
for caseIndex = 1:3
    detail = representativeDetails{caseIndex};
    points = detail.validationPoints{1};
    colour = log10(max(detail.validationRelative{1},eps));
    ax = nexttile(residualLayout,caseIndex);
    scatter3(ax,points(:,1),points(:,2),points(:,3),4,colour,'filled');
    axis(ax,'equal'); view(ax,3); colorbar(ax);
    title(ax,labels(caseIndex)+": full sphere");
    ax = nexttile(residualLayout,3+caseIndex);
    scatter3(ax,points(:,1),points(:,2),points(:,3),7,colour,'filled');
    axis(ax,'equal'); view(ax,3); colorbar(ax);
    xlim(ax,[0.9,1.01]); ylim(ax,[-0.16,0.16]);
    zlim(ax,[-0.16,0.16]);
    title(ax,labels(caseIndex)+": gap, log10 residual");

    collocation = detail.grid.body(1).collocation;
    caps = detail.grid.body(1).capPoints;
    ax = nexttile(distributionLayout,caseIndex);
    scatter3(ax,collocation(:,1),collocation(:,2), ...
        collocation(:,3),3,[0.5,0.5,0.5],'filled');
    hold(ax,'on');
    if ~isempty(caps)
        scatter3(ax,caps(:,1),caps(:,2),caps(:,3), ...
            8,[0.1,0.3,0.9],'filled');
    end
    axis(ax,'equal'); view(ax,3);
    title(ax,labels(caseIndex));

    ax = nexttile(closeLayout,caseIndex);
    facing = collocation(:,1) > 0.9;
    scatter(ax,collocation(facing,2),collocation(facing,3), ...
        3,[0.5,0.5,0.5],'filled');
    hold(ax,'on');
    if ~isempty(caps)
        scatter(ax,caps(:,2),caps(:,3),8,[0.1,0.3,0.9],'filled');
    end
    scatter(ax,0,0,40,'k','filled');
    axis(ax,'equal'); grid(ax,'on');
    xlim(ax,[-0.15,0.15]); ylim(ax,[-0.15,0.15]);
    xlabel(ax,'Tangent y'); ylabel(ax,'Tangent z');
    title(ax,labels(caseIndex)+": gap");
end
