function results = oct27_brownian_stress_divergence_images()
%OCT27_BROWNIAN_STRESS_DIVERGENCE Fixed-N, image-enriched two-sphere stress.
% Run deliberately: each gap needs 13 solves (25 with the half-step check).
% Unit spheres, unit viscosity, no imposed flow. 

%% Controls
N = 686;                         % Fixed coarse source count per sphere.
gapRatios = logspace(-2,0,9);     % Surface separation delta/a; a = 1.
imageFamilies = "SRD";            % S: Stokeslets, R: rotlets, D: source dipoles.
nearThreshold = 0.15;
proxyRadius = 0.63;
collocationRatio = 1.2;
thermalEnergy = 1;               % k_B T; ratio scales linearly with this.
volume = 1;                      % Common stress normalization; cancels in ratio.
forceMagnitude = 1;              % Each body is pulled toward its neighbour.
stepFraction = 0.02;             % h = min(0.01*a, stepFraction*delta).
checkHalfStep = true;            % 12 extra solves per gap when enabled.

root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(root,'setpath.m'));
opt = image_solver_options(struct('N',N,'proxyRadius',proxyRadius, ...
    'collocationRatio',collocationRatio,'imageFamilies',imageFamilies, ...
    'nearThreshold',nearThreshold,'gmresTolerance',1e-9, ...
    'maxIterations',400));

%% Fixed coarse grid, gap-dependent enrichment
% Preserve D_p = sum_(body j, axis alpha) d(S_p M e_jalpha)/d q_jalpha.
% A unit-force mobility solve evaluates the complete product S M directly.
% For spheres only translations contribute: physical orientation is immaterial.
% The paper rule counts spatial image nodes per body, shared by all families.
imageCounts = max(0,ceil(-8.72*log10(gapRatios)-6.15));
imageCounts(gapRatios >= nearThreshold) = 0;
nGaps = numel(gapRatios);
actualN = zeros(nGaps,1);
collocationCount = zeros(nGaps,1);
translationSteps = zeros(nGaps,1);
divergence = zeros(3,3,2,nGaps);
deterministicStresslet = zeros(3,3,2,nGaps);
deterministicStress = zeros(3,3,nGaps);
brownianStress = zeros(3,3,nGaps);
stepChange = NaN(nGaps,1);
surfaceResidual = zeros(nGaps,1);
iterations = zeros(nGaps,1);
forceMomentError = zeros(nGaps,1);

for k = 1:nGaps
    delta = gapRatios(k);
    q = [-1-delta/2,0,0;1+delta/2,0,0];
    h = min(0.01,stepFraction*delta);
    translationSteps(k) = h;
    % Freeze both count and activation across +/-h, including threshold gaps.
    opt.N_image = imageCounts(k);
    opt.nearThreshold = Inf;

    % Body layout: [Fx;Fy;Fz;Tx;Ty;Tz] for body 1, then body 2.
    neighbourDirection = (q(2,:)-q(1,:))/norm(q(2,:)-q(1,:));
    F = zeros(12,1);
    F(1:3) = forceMagnitude*neighbourDirection.';
    F(7:9) = -F(1:3);
    [~,it,~,err,details] = solve_mobility_images(q,F,opt);
    deterministicStresslet(:,:,:,k) = ...
        apply_KS_transpose(details.grid,details.lambda);
    actualN(k) = details.grid.N;
    collocationCount(k) = size(details.grid.body(1).collocation,1);
    forceMomentError(k) = norm(details.realisedForceTorque-F,inf);

    [D,divErr,divIt,divForceError] = image_stress_divergence(q,opt,h);
    if checkHalfStep
        [Dhalf,halfErr,halfIt,halfForceError] = ...
            image_stress_divergence(q,opt,h/2);
        stepChange(k) = norm(Dhalf(:)-D(:))/max(norm(Dhalf(:)),eps);
        D = Dhalf;              % Use the finer estimate in the comparison.
        divErr = max(divErr,halfErr);
        divIt = max(divIt,halfIt);
        divForceError = max(divForceError,halfForceError);
    end
    divergence(:,:,:,k) = D;
    surfaceResidual(k) = max(err,divErr);
    iterations(k) = max(it,divIt);
    forceMomentError(k) = max(forceMomentError(k),divForceError);
    % Force-on-fluid moments: suspension stress carries a minus sign.
    % Compare hydrodynamic stresslets; no centre-force virial is included.
    deterministicStress(:,:,k) = -sum(deterministicStresslet(:,:,:,k),3)/volume;
    brownianStress(:,:,k) = -thermalEnergy*sum(D,3)/volume;
    fprintf(['delta/a %.4g, N %d, images/body %d: surface %.2e, ', ...
        'step change %.2e, max GMRES %d\n'],delta,actualN(k), ...
        imageCounts(k),surfaceResidual(k),stepChange(k),iterations(k));
end

%% Signed stress components and relative size
deterministicDiagonal = zeros(nGaps,3);
brownianDiagonal = zeros(nGaps,3);
relativeSize = zeros(nGaps,1);
for k = 1:nGaps
    deterministicDiagonal(k,:) = diag(deterministicStress(:,:,k)).';
    brownianDiagonal(k,:) = diag(brownianStress(:,:,k)).';
    relativeSize(k) = norm(brownianStress(:,:,k),'fro')/ ...
        max(norm(deterministicStress(:,:,k),'fro'),eps);
end
results = table(gapRatios(:),actualN,imageCounts(:),collocationCount, ...
    relativeSize,stepChange,surfaceResidual,iterations,forceMomentError, ...
    'VariableNames',{'deltaOverA','N','imagesPerBody','collocationPerBody', ...
    'brownianOverDeterministic','stepChange','surfaceResidual', ...
    'maxIterations','forceMomentError'});
disp(results);

fig = figure('Name','Oct 27: deterministic and Brownian stress', ...
    'Position',[100,100,1200,650]);
layout = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
componentNames = {'xx','yy','zz'};
for axis = 1:3
    ax = nexttile(layout);
    semilogx(ax,gapRatios,deterministicDiagonal(:,axis),'o-', ...
        'LineWidth',1.4,'DisplayName','Deterministic');
    hold(ax,'on');
    semilogx(ax,gapRatios,brownianDiagonal(:,axis),'s--', ...
        'LineWidth',1.4,'DisplayName','Brownian divergence');
    yline(ax,0,':','HandleVisibility','off');
    grid(ax,'on'); xlabel(ax,'Surface gap \delta/a');
    ylabel(ax,['\sigma_{',componentNames{axis},'}']);
    title(ax,['Signed ',componentNames{axis},' stress']);
    legend(ax,'Location','best');
end
ax = nexttile(layout);
loglog(ax,gapRatios,relativeSize,'o-','LineWidth',1.4);
hold(ax,'on'); yline(ax,1,':','Equal magnitude'); grid(ax,'on');
xlabel(ax,'Surface gap \delta/a');
ylabel(ax,'||\sigma_B||_F / ||\sigma_{det}||_F');
title(ax,sprintf('Relative size, k_B T = %g',thermalEnergy));
ax = nexttile(layout);
semilogx(ax,gapRatios,imageCounts,'o-','LineWidth',1.4);
grid(ax,'on'); xlabel(ax,'Surface gap \delta/a');
ylabel(ax,'Image nodes per body');
title(ax,sprintf('Fixed N = %d, images %s',actualN(1),imageFamilies));
ax = nexttile(layout);
loglog(ax,gapRatios,surfaceResidual,'o-','DisplayName','Surface residual');
hold(ax,'on');
if checkHalfStep
    loglog(ax,gapRatios,stepChange,'s--','DisplayName','h to h/2 change');
end
grid(ax,'on'); xlabel(ax,'Surface gap \delta/a'); ylabel(ax,'Relative error');
title(ax,'Numerical diagnostics'); legend(ax,'Location','best');
title(layout,'Two unit spheres: inward unit forces; hydrodynamic stresslets');

outputFolder = fullfile(root,'visuals');
if ~isfolder(outputFolder), mkdir(outputFolder); end
exportgraphics(fig,fullfile(outputFolder,'oct27_brownian_stress_comparison.png'));
save(fullfile(root,'experiments','oct27_brownian_stress_divergence_results.mat'), ...
    'results','divergence','deterministicStresslet','deterministicStress', ...
    'brownianStress','translationSteps','checkHalfStep','N','gapRatios', ...
    'imageCounts','opt','thermalEnergy','volume','forceMagnitude');
end
