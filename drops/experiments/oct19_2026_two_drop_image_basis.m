% OCT19_2026_TWO_DROP_IMAGE_BASIS Separate direction and outer-line effects.
% A deliberately small dense reference: geometry, forcing, collocation and
% validation nodes stay fixed. Dense gap caps and a cap-refinement check
% distinguish cap resolution from basis effects. Increase N separately to
% test a finer proxy grid; the default is a small source discretisation.
clear; close all;
plot_grid = true;
runIterative = true;
capCounts = [1200 2400];       % First run, then the cap-resolution check.
capLayers = [40 57];          % Refine both meridional and azimuthal sampling.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(root));
q = [-1.05 0 0;1.05 0 0]; a = ones(2,1);
S1 = ones(2,1); e = repmat([0 0 1],2,1);
% imageFamilies tokens: S Stokeslet, R rotlet, D potential dipole, T legacy
% stresslet, T_t line-parallel, T_n first transverse, T_n2 second transverse
% stresslet, Q_s trace-free axial Stokeslet quadrupole. T excludes T_*.
opt = struct('N',100,'collocationRatio',2,'N_image',4, ...
    'imageFamilies',"SDRT_nT_t",'mirror_ext',true, ...
    'nearThreshold',0.15,'proxyRadius',0.63,'outerProxyRatio',2, ...
    'outerImageOffsetFraction',0.01,'capMode',"exponential", ...
    'capPointCounts',capCounts(1),'capLayers',capLayers(1), ...
    'useRowScaling',false,'interiorViscosity',2, ...
    'svdCutoff',1e-10,'gmresTolerance',1e-10,'maxIterations',150, ...
    'validationCount',2000,'validationGapCount',4001,'plot_grid',plot_grid, ...
    'solverMode',"denseSVD");
opt = drop_image_options(opt,2,a);
base = drop_image_discretization(q,a,opt);
[points,normals] = drop_image_validation_grid(base, ...
    opt.validationCount,opt.validationGapCount,pi/5);
collocation = {base.body.collocation};

%% First vary the basis, then exterior-image presence and endpoint only
mixes = ["SDR","SDRT_t","SDRT_n","SDRT_nT_t", ...
    "SDRT_nT_t","SDRT_nT_t","SDRT_nT_t"];
mirrors = [true true true true false true true];
offsets = [0.01 0.01 0.01 0.01 0.01 0.2 0.005];
labels = ["SDR","SDRT_t","SDRT_n","SDRT_nT_t", ...
    "joint inner only","joint endpoint 1.60a","joint endpoint 1.99a"];
families = image_family_list("SDRT_nT_t");
nr = numel(labels);
residual = zeros(nr,3); gapError = residual; awayError = residual;
proxy = zeros(nr,2); innerImages = nan(nr,5); outerImages = innerImages;
rankKept = zeros(nr,1); algebraic = rankKept; cols = rankKept;
translations = zeros(2,3,nr);
pointErrors = cell(nr,1);
scale = [max(abs(S1).*a),max(abs(S1).*a),max(abs(S1))];
for k = 1:nr
    opt.imageFamilies = mixes(k); opt.mirror_ext = mirrors(k); % See token legend above.
    opt.outerImageOffsetFraction = offsets(k);
    fprintf('\n%s\n',labels(k));
    geometry = drop_image_discretization(q,a,opt);
    disp(drop_image_check_geometry(geometry,opt));
    [U,lambda,gamma,g,info] = solve_drops_with_images(q,a,S1,e,opt);
    for i = 1:2
        if ~isequal(g.body(i).collocation,collocation{i})
            error('oct19_2026_two_drop_image_basis:ChangedGrid','Collocation changed.');
        end
    end
    [maxError,pointError] = drop_image_validate_surface_conditions( ...
        g,lambda,gamma,U,S1,e,info.options,points,normals);
    residual(k,:) = max(maxError,[],1)./scale;
    for i = 1:2
        ray = g.body(i).contactDirections(1,:);
        inGap = normals{i}*ray.'>=cos(pi/5);
        gapError(k,:) = max(gapError(k,:),max(pointError{i}(inGap,:),[],1)./scale);
        awayError(k,:) = max(awayError(k,:),max(pointError{i}(~inGap,:),[],1)./scale);
    end
    [proxy(k,1),proxy(k,2),innerImages(k,:),outerImages(k,:)] = ...
        drop_image_max_source_strengths(g,lambda,gamma,families);
    rankKept(k) = info.retainedRanks; algebraic(k) = info.originalResidual;
    cols(k) = info.systemSize(2); translations(:,:,k) = U;
    pointErrors{k} = pointError;
    if k==4
        drop_image_plot_source_strengths(g,lambda,gamma,'SDRT_nT_t dense');
    end
end

%% Tabulate independent errors, truncation, and strengths by source family
results = table(labels.',residual(:,1),residual(:,2),residual(:,3), ...
    rankKept,cols,algebraic,'VariableNames',{'Case','Continuity', ...
    'NormalVelocity','TangentialTraction','Rank','Columns','AlgebraicResidual'});
disp(results);
regionNames = {'Continuity','NormalVelocity','TangentialTraction'};
fprintf('Gap-facing cap maxima (angular radius pi/5):\n');
disp(array2table(gapError,'VariableNames',regionNames,'RowNames',cellstr(labels)));
fprintf('Away-from-gap maxima:\n');
disp(array2table(awayError,'VariableNames',regionNames,'RowNames',cellstr(labels)));
fprintf('Inner source magnitudes: proxy then each image family (NaN absent):\n');
disp(array2table([proxy(:,1),innerImages], ...
    'VariableNames',cellstr(["ProxyS",families]),'RowNames',cellstr(labels)));
fprintf('Outer source magnitudes: proxy then each image family:\n');
disp(array2table([proxy(:,2),outerImages], ...
    'VariableNames',cellstr(["ProxyS",families]),'RowNames',cellstr(labels)));

%% Compare the iterative method on exactly the fourth dense discretisation
if runIterative
    opt.imageFamilies = mixes(4); opt.mirror_ext = mirrors(4); % See token legend above.
    opt.outerImageOffsetFraction = offsets(4); opt.solverMode = "oneBody";
    [U,lambda,gamma,g,info] = solve_drops_with_images(q,a,S1,e,opt);
    maxError = drop_image_validate_surface_conditions( ...
        g,lambda,gamma,U,S1,e,info.options,points,normals);
    fprintf(['Joint iterative: iterations %d, original residual %.3e, ', ...
        'preconditioned residual %.3e\n'], ...
        info.iterations,info.originalResidual,info.regularizedResidual);
    fprintf('Joint iterative independent maxima: %.3e %.3e %.3e\n', ...
        max(maxError,[],1)./scale);
    fprintf('Dense/iterative absolute translation difference: %.3e\n', ...
        norm(U-translations(:,:,4),'fro'));
end

%% Double the dense cap for SDR and the joint basis on the same validation grid
refinedResidual = zeros(2,3); refinedGap = refinedResidual; refinedAway = refinedResidual;
refinedAlgebraic = zeros(2,1); refinedRanks = zeros(2,1);
for k = 1:2
    index = [1 4]; index = index(k);
    checkOpt = opt;
    checkOpt.solverMode = "denseSVD";
    checkOpt.imageFamilies = mixes(index); checkOpt.mirror_ext = true; % See token legend above.
    checkOpt.outerImageOffsetFraction = 0.01;
    checkOpt.capPointCounts = capCounts(2); checkOpt.capLayers = capLayers(2);
    fprintf('\nCap refinement: %s, requested cap count %d\n',mixes(index),capCounts(2));
    [U,lambda,gamma,g,info] = solve_drops_with_images(q,a,S1,e,checkOpt);
    [maxError,pointError] = drop_image_validate_surface_conditions( ...
        g,lambda,gamma,U,S1,e,info.options,points,normals);
    refinedResidual(k,:) = max(maxError,[],1)./scale;
    for i = 1:2
        inGap = normals{i}*g.body(i).contactDirections(1,:).'>=cos(pi/5);
        refinedGap(k,:) = max(refinedGap(k,:),max(pointError{i}(inGap,:),[],1)./scale);
        refinedAway(k,:) = max(refinedAway(k,:),max(pointError{i}(~inGap,:),[],1)./scale);
    end
    refinedAlgebraic(k) = info.originalResidual; refinedRanks(k) = info.retainedRanks;
    fprintf('Actual cap count on drop 1: %d\n',size(g.body(1).capPoints,1));
end
fprintf('Doubled-cap independent maxima:\n');
disp(array2table(refinedResidual,'VariableNames',regionNames, ...
    'RowNames',cellstr(mixes([1 4]))));
fprintf('Doubled-cap gap maxima:\n');
disp(array2table(refinedGap,'VariableNames',regionNames, ...
    'RowNames',cellstr(mixes([1 4]))));
fprintf('Doubled-cap away-from-gap maxima:\n');
disp(array2table(refinedAway,'VariableNames',regionNames, ...
    'RowNames',cellstr(mixes([1 4]))));
% Rows remain unweighted, so changing cap density also changes its relative
% least-squares weight. Compare independent residuals, not algebraic norms alone.

%% Plot the three conditions globally and split by gap location
figure('Name','Drop basis and outer-line comparison','Color','w', ...
    'Position',[100 100 1500 540]);
tiledlayout(1,3,'TileSpacing','compact');
conditionNames = {'Continuity','Normal velocity','Tangential traction'};
for j = 1:3
    ax = nexttile;
    semilogy(ax,1:nr,residual(:,j),'-o',1:nr,gapError(:,j),'--s', ...
        1:nr,awayError(:,j),':^','LineWidth',1.2);
    hold(ax,'on');
    semilogy(ax,[1 4],refinedResidual(:,j),'kx','MarkerSize',9,'LineWidth',1.6);
    xticks(ax,1:nr); xticklabels(ax,labels); xtickangle(ax,35);
    set(ax,'XGrid','on','YGrid','on');
    title(ax,conditionNames{j}); ylabel(ax,'Maximum scaled residual');
end
lg = legend(ax,{'whole surface','gap cap','away from gap','whole, doubled cap'}, ...
    'Orientation','horizontal');
lg.Layout.Tile = 'south';

% A common colour range makes residual location comparable across bases.
figure('Name','Surface residuals for SDR and the joint stresslets','Color','w', ...
    'Position',[100 100 1400 850]);
tiledlayout(2,3,'TileSpacing','compact');
for k = [1 4]
    for j = 1:3
        ax = nexttile; hold(ax,'on');
        for i = 1:2
            value = log10(max(pointErrors{k}{i}(:,j)/scale(j),realmin));
            scatter3(ax,points{i}(:,1),points{i}(:,2),points{i}(:,3),9,value,'filled');
        end
        axis(ax,'equal'); view(ax,3); clim(ax,[-6 0]); colorbar(ax);
        title(ax,labels(k)+": "+regionNames{j});
    end
end
colormap(turbo);
