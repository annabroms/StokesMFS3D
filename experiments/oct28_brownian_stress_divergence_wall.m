%OCT28_BROWNIAN_STRESS_DIVERGENCE_WALL One force-driven sphere above a wall.
% Deliberate execution: 3 solves/height, or 5 with the half-step check.
% Preserve the complete unit-force stresslet response when differentiating.

%% Controls
radius = 1;
viscosity = 1;
wallHeight = 0;
gapRatios = logspace(-3,0,13);      % delta/a; centre height above wall h=a+delta.
N = 700;                          % Spherical designs may round down.
proxyRadius = 0.4;                % Proxy radius / a.
collocationRatio = 3;
wallImageFamilies = "ST_tD";
imageCount = 30;
imageGapThreshold = 0.2;           % Strict bound on delta/a.
imageOffsetFraction = 0.05;
capAngle = pi/5;
capCount = 901;
capLayers = 30;
capExponent = 4;
svdCutoff = 1e-12;
SVDsolve = false;
useFmm = false;
gmresTolerance = 1e-9;
fmmTolerance = 1e-12;
maxIterations = 400;
validationResolution = [32,40];
validationCount = 2001;
batchSize = 256;
thermalEnergy = 1;                % k_B T, physical energy units.
volume = 1;                       % Common normalization; cancels in ratio.
forceMagnitude = 1;               % Applied force [0;0;-forceMagnitude].
stepFraction = 0.02;              % Initial epsilon <= stepFraction*delta.
maximumStepRatio = 0.01;          % Initial epsilon <= maximumStepRatio*a.
checkHalfStep = true;
stressletCheckResolution = [48,64]; % [] disables base-solve reciprocity diagnostic.
plotAgainstGap = false;           % false: centre h/a; true: clearance delta/a.
outputPrefix = 'oct28_brownian_stress_wall';

%% Geometry limits and complete solver options
projectRoot = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(projectRoot,'setpath.m'));
validateattributes(gapRatios,{'double'},{'vector','real','finite','positive'});
validateattributes(radius,{'double'},{'scalar','real','finite','positive'});
validateattributes(volume,{'double'},{'scalar','real','finite','positive'});
validateattributes(thermalEnergy,{'double'},{'scalar','real','finite','nonnegative'});
validateattributes(forceMagnitude,{'double'},{'scalar','real','finite','positive'});
validateattributes(stepFraction,{'double'},{'scalar','real','finite','>',0,'<',1});
validateattributes(maximumStepRatio,{'double'},{'scalar','real','finite','positive'});
validateattributes(imageGapThreshold,{'double'},{'scalar','real','finite','positive'});
rho = proxyRadius*(1+imageOffsetFraction);
validateattributes(rho,{'double'},{'scalar','real','finite','>',0,'<',1});
% rho < h/a-sqrt((h/a)^2-1) iff delta/a < (1-rho)^2/(2*rho).
geometricMaxGapRatio = (1-rho)^2/(2*rho);
maxEnhancedGapRatio = min(imageGapThreshold,geometricMaxGapRatio);
fprintf('Wall images require delta/a < %.8g; geometric limit %.8g.\n', ...
    maxEnhancedGapRatio,geometricMaxGapRatio);
solverOptions = struct('radius',radius,'viscosity',viscosity, ...
    'wall_height',wallHeight,'wallImageFamilies',wallImageFamilies, ...
    'pairImageFamilies',"none",'N_image',imageCount, ...
    'imageOffsetFraction',imageOffsetFraction,'capAngle',capAngle, ...
    'capCount',capCount,'capLayers',capLayers,'capExponent',capExponent, ...
    'tol',svdCutoff,'SVDsolve',SVDsolve,'fmm',useFmm, ...
    'gmres_tol',gmresTolerance,'fmm_tol',fmmTolerance,'maxit',maxIterations, ...
    'gmres_verbose',0,'compute_residual',1,'plot',0, ...
    'validationResolution',validationResolution,'validationCount',validationCount, ...
    'batchSize',batchSize);
discretization = struct('N',N,'proxyRadius',proxyRadius, ...
    'collocationRatio',collocationRatio);
parameters = struct('radius',radius,'viscosity',viscosity,'wallHeight',wallHeight, ...
    'gapRatios',gapRatios,'N',N,'proxyRadius',proxyRadius, ...
    'collocationRatio',collocationRatio,'wallImageFamilies',wallImageFamilies, ...
    'imageCount',imageCount,'imageGapThreshold',imageGapThreshold, ...
    'imageOffsetFraction',imageOffsetFraction,'thermalEnergy',thermalEnergy, ...
    'volume',volume,'forceMagnitude',forceMagnitude,'stepFraction',stepFraction, ...
    'maximumStepRatio',maximumStepRatio,'checkHalfStep',checkHalfStep, ...
    'stressletCheckResolution',stressletCheckResolution,'plotAgainstGap',plotAgainstGap);

gapRatios = gapRatios(:); nHeights = numel(gapRatios);
centreHeights = radius*(1+gapRatios);
divergence = zeros(3,3,nHeights);
deterministicStresslet = zeros(3,3,nHeights);
deterministicStress = zeros(3,3,nHeights);
brownianStress = zeros(3,3,nHeights);
translationSteps = zeros(nHeights,1);
diagnostics = cell(nHeights,1);
relativeSize = NaN(nHeights,1);
summary = zeros(nHeights,11);

%% Solve and differentiate with a fixed image basis across each stencil
for k = 1:nHeights
    opt = solverOptions;
    if gapRatios(k) >= maxEnhancedGapRatio
        opt.N_image = 0;
    end
    step = radius*min(maximumStepRatio,stepFraction*gapRatios(k));
    if opt.N_image > 0 && ~strcmpi(wallImageFamilies,"none")
        % Stay strictly inside the geometric bound; never differentiate an
        % automatic source-count switch. Activation at 0.2 is frozen locally.
        step = min(step,0.25*radius*(geometricMaxGapRatio-gapRatios(k)));
    end
    if step <= 100*eps(max(abs(wallHeight)+centreHeights(k),radius))
        error('wall:StressStep','Derivative step is unresolvable; move the sample off the cutoff.');
    end
    translationSteps(k) = step;
    q = [0,0,wallHeight+centreHeights(k)];
    [D,Sdown,info] = single_sphere_wall_stress_divergence( ...
        q,discretization,opt,step,checkHalfStep,stressletCheckResolution);
    divergence(:,:,k) = D;
    deterministicStresslet(:,:,k) = forceMagnitude*Sdown;
    % Force-on-fluid stresslets have the same suspension-stress minus sign
    % as the Oct 27 experiment. No centre-force virial or isotropic pressure.
    deterministicStress(:,:,k) = -forceMagnitude*Sdown/volume;
    brownianStress(:,:,k) = -thermalEnergy*D/volume;
    denominator = norm(deterministicStress(:,:,k),'fro');
    if denominator > 0
        relativeSize(k) = norm(brownianStress(:,:,k),'fro')/denominator;
    end
    diagnostics{k} = info;
    summary(k,:) = [info.proxyCount,info.imagePoints,info.sourceVectors, ...
        info.collocationPoints,max(info.surfaceResidual),max(info.momentError), ...
        max(info.iterations),info.stepChange,info.stressletCorrection, ...
        max(info.rank)-min(info.rank),max(info.linearSystemResidual)];
    fprintf(['h/a %.6g, delta/a %.4g, images %d: surface %.3e, ', ...
        'step change %.3e, Brownian/deterministic %.3e\n'], ...
        centreHeights(k)/radius,gapRatios(k),info.imagePoints, ...
        summary(k,5),info.stepChange,relativeSize(k));
end
results = array2table(summary,'VariableNames',{'proxyCount','imagePoints', ...
    'sourceVectors','collocationPoints','surfaceResidual','momentError', ...
    'maxIterations','stepChange','stressletCorrection','rankVariation', ...
    'linearSystemResidual'});
results = addvars(results,centreHeights/radius,gapRatios,translationSteps,relativeSize, ...
    'Before',1,'NewVariableNames',{'heightOverA','gapOverA','step','brownianOverDeterministic'});
divergenceDiagonal = zeros(nHeights,3);
deterministicDiagonal = zeros(nHeights,3);
brownianDiagonal = zeros(nHeights,3);
for k = 1:nHeights
    divergenceDiagonal(k,:) = diag(divergence(:,:,k)).';
    deterministicDiagonal(k,:) = diag(deterministicStress(:,:,k)).';
    brownianDiagonal(k,:) = diag(brownianStress(:,:,k)).';
end
% Save signed components in both MAT and CSV, even though figures show magnitudes.
components = {'xx','yy','zz'};
for j = 1:3
    results.(['D_',components{j}]) = divergenceDiagonal(:,j);
    results.(['sigmaDet_',components{j}]) = deterministicDiagonal(:,j);
    results.(['sigmaB_',components{j}]) = brownianDiagonal(:,j);
end
disp(results);

%% Log-log visualization; crosses identify negative signed components
abscissa = centreHeights/radius;
axisLabel = 'Sphere centre height above wall, h/a';
cutoff = 1+maxEnhancedGapRatio;
if plotAgainstGap
    abscissa = gapRatios; %#ok<UNRCH> Selected by the top-level control.
    axisLabel = 'Sphere-wall surface gap, \delta/a';
    cutoff = maxEnhancedGapRatio;
end
fig = figure('Name','Oct 28: wall Brownian stress divergence', ...
    'Position',[100,100,1250,850]);
layout = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
colors = lines(3);
ax = nexttile(layout); hold(ax,'on');
for j = 1:3
    signed_magnitude(ax,abscissa,divergenceDiagonal(:,j),colors(j,:), ...
        '-o',['D_{',components{j},'}']);
end
ylabel(ax,'|D_{ii}|'); title(ax,'Stresslet-response divergence');
legend(ax,'Location','best');
ax = nexttile(layout); hold(ax,'on');
for j = 1:3
    signed_magnitude(ax,abscissa,deterministicDiagonal(:,j),colors(j,:), ...
        '-o',['Det ',components{j}]);
    signed_magnitude(ax,abscissa,brownianDiagonal(:,j),colors(j,:), ...
        '--s',['Brownian ',components{j}]);
end
ylabel(ax,'|\sigma_{ii}|'); title(ax,'Deterministic and Brownian stresses');
legend(ax,'Location','best');
ax = nexttile(layout);
loglog(ax,abscissa,relativeSize,'o-','LineWidth',1.4);
yline(ax,1,':','Equal magnitude');
ylabel(ax,'||\sigma_B||_F / ||\sigma_{det}||_F'); title(ax,'Relative stress magnitude');
ax = nexttile(layout);
loglog(ax,abscissa,results.surfaceResidual,'o-','DisplayName','Surface residual');
hold(ax,'on');
if checkHalfStep
    loglog(ax,abscissa,results.stepChange,'s--','DisplayName','Step-halving change');
end
loglog(ax,abscissa,results.stressletCorrection,'d:', ...
    'DisplayName','Reciprocal stresslet correction');
ylabel(ax,'Relative diagnostic'); title(ax,'Accuracy diagnostics');
legend(ax,'Location','best');
axesList = findall(fig,'Type','axes');
for j = 1:numel(axesList)
    ax = axesList(j);
    set(ax,'XScale','log','YScale','log'); grid(ax,'on'); xlabel(ax,axisLabel);
    xline(ax,cutoff,':','HandleVisibility','off');
end
title(layout,sprintf(['Single sphere: F_z=-%g, k_BT=%g, images %s; ', ...
    'crosses mark negative components, zeros omitted'], ...
    forceMagnitude,thermalEnergy,wallImageFamilies),'Interpreter','none');

%% Save full signed data and reproducible settings
outputDirectory = fullfile(projectRoot,'visuals');
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
outputStem = fullfile(outputDirectory,outputPrefix);
save([outputStem,'.mat'],'results','divergence','deterministicStresslet', ...
    'deterministicStress','brownianStress','translationSteps','diagnostics', ...
    'parameters','solverOptions','geometricMaxGapRatio','maxEnhancedGapRatio');
writetable(results,[outputStem,'.csv']);
savefig(fig,[outputStem,'.fig']);
exportgraphics(fig,[outputStem,'.png'],'Resolution',200);
exportgraphics(fig,[outputStem,'.pdf'],'ContentType','vector');

function signed_magnitude(ax,x,y,color,style,label)
% Keep log-log magnitudes readable while marking negative values explicitly.
magnitude = abs(y); magnitude(magnitude==0) = NaN;
loglog(ax,x,magnitude,style,'Color',color,'LineWidth',1.3,'DisplayName',label);
negative = y<0;
loglog(ax,x(negative),magnitude(negative),'x','Color',color, ...
    'MarkerSize',9,'LineWidth',1.5,'HandleVisibility','off');
end
