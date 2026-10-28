%OCT28_SINGLE_SPHERE_WALL_MOBILITY Unit downward force at varying wall heights.
% Deliberate experiment: one image-enriched mobility solve per height.
% Preserve F=[0;0;-1], zero torque and the complete Blake wall correction
% throughout the sweep; report the independently sampled surface residual.

%% Controls
radius = 1;
viscosity = 1;
wallHeight = 0;
gapRatios = logspace(-1,0,5);     % Surface gap delta/a; centre height h=a+delta.
forceTorque = [0;0;-1;0;0;0];     % Physical applied load; solver's on-fluid sign.
N = 700;                         % Requested proxy count; design may round down.
proxyRadius = 0.4;               % Proxy radius divided by sphere radius.
collocationRatio = 3;             % Keeps proxy-only fits >=2 rows per unknown.
wallImageFamilies = "ST_tD";
imageCount = 30;                  % Geometric nodes per active wall line.
imageGapThreshold = 0.2;          % Enable enrichment only for delta/a < this.
imageOffsetFraction = 0.05;
capAngle = pi/5;
capCount = 901;
capLayers = 30;
capExponent = 4;
svdCutoff = 1e-12;
SVDsolve = false;                 % true: global rectangular SVD reference.
useFmm = false;                   % true requires compiled FMM3D for S sources.
gmresTolerance = 1e-9;
fmmTolerance = 1e-12;
maxIterations = 400;
validationResolution = [32,40];
validationCount = 2001;           % Per independent wall-facing validation cap.
batchSize = 256;
plotAgainstGap = false;           % false: h/a; true: surface clearance delta/a.

%% Largest wall gap eligible for image enhancement
% Rp = proxyRadius is dimensionless (physical proxy radius divided by a).
% Let rho = Rp*(1+imageOffsetFraction), the inner image-line cutoff / a.
% At centre height H=h/a, the accumulation radius / a is
% rstar = H-sqrt(H^2-1). A nonempty line requires rho < rstar.
% For 0 < rho < 1 this is equivalent to
% delta/a < (rho+1/rho)/2-1 = (1-rho)^2/(2*rho).
% With zero offset: delta_max/a = (1-Rp)^2/(2*Rp).
% Both this geometric bound and imageGapThreshold are STRICT upper bounds;
% equality has no enhancement. Positive imageCount and active families are
% also required. These are sphere-wall gaps, not reflected-pair gaps (2*delta).
imageInnerRadiusRatio = proxyRadius*(1+imageOffsetFraction);
validateattributes(imageInnerRadiusRatio,{'double'}, ...
    {'scalar','real','finite','>',0,'<',1});
geometricMaxGapRatio = (1-imageInnerRadiusRatio)^2/(2*imageInnerRadiusRatio);
maxEnhancedGapRatio = min(imageGapThreshold,geometricMaxGapRatio);
maxEnhancedGap = radius*maxEnhancedGapRatio;
maxEnhancedCentreHeight = radius+maxEnhancedGap; % Relative to wall plane.
fprintf(['Image enhancement requires delta/a < %.8g ', ...
    '(geometry %.8g, configured threshold %.8g); h/a < %.8g.\n'], ...
    maxEnhancedGapRatio,geometricMaxGapRatio,imageGapThreshold, ...
    maxEnhancedCentreHeight/radius);

%% Setup
projectRoot = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(projectRoot,'setpath.m'));
validateattributes(gapRatios,{'double'},{'vector','real','finite','positive'});
validateattributes(imageGapThreshold,{'double'},{'scalar','real','finite','positive'});
opt = struct('radius',radius,'viscosity',viscosity,'wall_height',wallHeight, ...
    'wallImageFamilies',wallImageFamilies,'pairImageFamilies',"none", ...
    'N_image',imageCount,'imageOffsetFraction',imageOffsetFraction, ...
    'capAngle',capAngle,'capCount',capCount,'capLayers',capLayers, ...
    'capExponent',capExponent,'tol',svdCutoff,'SVDsolve',SVDsolve, ...
    'fmm',useFmm,'gmres_tol',gmresTolerance,'fmm_tol',fmmTolerance, ...
    'maxit',maxIterations,'gmres_verbose',0,'compute_residual',1,'plot',0, ...
    'validationResolution',validationResolution,'validationCount',validationCount, ...
    'batchSize',batchSize);

gapRatios = gapRatios(:);
wallGaps = radius*gapRatios;
centreHeights = radius+wallGaps;    % Measured relative to the wall plane.
nHeights = numel(gapRatios);
rigidVelocity = zeros(6,nHeights);
surfaceResidual = zeros(nHeights,1);
collocationResidual = zeros(nHeights,1);
iterations = zeros(nHeights,1);
momentError = zeros(nHeights,1);
linearSystemResidual = zeros(nHeights,1);
actualN = zeros(nHeights,1);
actualImageCount = zeros(nHeights,1);
sourceCount = zeros(nHeights,1);
collocationCount = zeros(nHeights,1);

%% Height sweep -- rebuild all wall-dependent geometry for each height
for k = 1:nHeights
    q = [0,0,wallHeight+centreHeights(k)];
    [rin,rout] = init_spheres(q,proxyRadius,N,collocationRatio,radius);
    opt.N_image = imageCount;
    if gapRatios(k) >= imageGapThreshold
        opt.N_image = 0;
    end
    [rigidVelocity(:,k),iterations(k),~,surfaceResidual(k),~,details] = ...
        solve_mobility_wall_images(q,rin,rout,forceTorque,opt);
    actualN(k) = details.N;
    sourceCount(k) = details.N_total;
    collocationCount(k) = details.M;
    collocationResidual(k) = details.collocationError;
    momentError(k) = details.momentError;
    linearSystemResidual(k) = details.linearSystemResidual;
    % The solver also omits lines whose accumulation point is inside the
    % proxy cutoff. Count actual geometry rather than the requested nodes.
    basis = details.problem.body{1}.basis;
    for b = 1:numel(basis)
        if startsWith(basis(b).label,"wall ")
            actualImageCount(k) = size(basis(b).points,1);
            break
        end
    end
    fprintf(['h/a %.6g, delta/a %.4g, images %d: ', ...
        'surface %.3e, moment error %.3e, iterations %d\n'], ...
        centreHeights(k)/radius,gapRatios(k),actualImageCount(k), ...
        surfaceResidual(k),momentError(k),iterations(k));
end

results = table(centreHeights/radius,gapRatios,actualN,actualImageCount, ...
    sourceCount,collocationCount,rigidVelocity(3,:).',surfaceResidual, ...
    collocationResidual,momentError,iterations,linearSystemResidual, ...
    'VariableNames',{'heightOverA','gapOverA','proxyCount','imagePoints', ...
    'sourceVectors','collocationPoints','verticalVelocity','surfaceResidual', ...
    'collocationResidual','momentError','iterations','linearSystemResidual'});
disp(results);
clear details                         % Release last geometry's SVD factors.

%% Independent relative surface residual versus height
fig = figure('Name','Oct 28: single-sphere wall mobility');
ax = axes(fig);
if plotAgainstGap
    abscissa = gapRatios; % Enabled by the top-level plotting control.
    axisLabel = 'Sphere-wall surface gap, \delta/a';
else
    abscissa = centreHeights/radius;
    axisLabel = 'Sphere centre height above wall, h/a';
end
loglog(ax,abscissa,surfaceResidual,'o-','LineWidth',1.4);
set(ax,'XScale','log','YScale','log');
grid(ax,'on');
xlabel(ax,axisLabel);
ylabel(ax,'Maximum relative surface velocity residual');
title(ax,sprintf('Oct 28: F_z = %g, wall images %s', ...
    forceTorque(3),wallImageFamilies));
% Residual = max independent |u-U-Omega cross (x-q)| / max sampled rigid speed.
% No plotting floor is applied; results retains the measured residuals.
