% OCT18_2026_TWO_DROP_IMAGE_COUNT_SWEEP Fixed two-drop image-resolution study.
% Run this script from any directory. Only N_image changes between solves.
% Counts start at one because zero images remove the cap collocation grid.
% Strength maxima are Euclidean norms of the physical three-component
% coefficients at one source point, maximised over both drops.

%% Fixed geometry, physics, discretisation, and solver controls
close all; clear;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(root));
q = [-1.05 0 0;1.05 0 0];
a = ones(2,1);
S1 = ones(2,1);
e = repmat([0 0 1],2,1);
imageCounts = [0 10];
% imageFamilies tokens: S Stokeslet, R rotlet, D potential dipole, T legacy
% stresslet, T_t line-parallel, T_n first transverse, T_n2 second transverse
% stresslet, Q_s trace-free axial Stokeslet quadrupole. T excludes T_*.
opt = struct('N',700,'collocationRatio',2,'nearThreshold',0.15, ...
    'N_image',imageCounts(1),'imageFamilies',"SRDT_tT_nT_n2Q_s",'useImages',true, ...
    'use_orth_dir_images',false, ...
    'mirror_ext',true,'proxyRadius',0.63,'outerProxyRatio',1.8, ...
    'outerImageOffsetFraction',0.01,'interiorViscosity',2, ...
    'svdCutoff',1e-12,'gmresTolerance',1e-9, ...
    'maxIterations',350,'batchSize',128,'solverMode',"oneBody", ...
    'validationCount',1000,'validationGapCount',401, ...
    'capMode',"exponential",'capPointCounts',1000, ...
    'useRowScaling',true,'plot_grid',false);
families = image_family_list(opt.imageFamilies);
nCases = numel(imageCounts);
nFamilies = numel(families);

% Generate one independent, gap-resolving validation grid for every count.
baseGrid = drop_image_discretization(q,a,drop_image_options(opt,2,a));
collocation = {baseGrid.body.collocation};
[points,normals] = drop_image_validation_grid( ...
    baseGrid,opt.validationCount,opt.validationGapCount,pi/5);
scale = [max(abs(S1).*a),max(abs(S1).*a),max(abs(S1))];
if any(scale==0)
    error('oct18_2026_two_drop_image_count_sweep:ZeroScale', ...
        'The active forcing must give positive residual scales.');
end

%% Solve each image resolution and measure independent surface errors
residual = zeros(nCases,3);
innerProxy = zeros(nCases,1);
outerProxy = zeros(nCases,1);
innerImages = nan(nCases,nFamilies);
outerImages = nan(nCases,nFamilies);
innerAll = zeros(nCases,1);
outerAll = zeros(nCases,1);
iterations = zeros(nCases,1);
algebraicResidual = zeros(nCases,1);
for k = 1:nCases
    opt.N_image = imageCounts(k);
    fprintf('\nImage nodes per line: %d\n',opt.N_image);
    [U,lambda,gamma,grid,info] = solve_drops_with_images( ...
        q,a,S1,e,opt);
    % for body = 1:2
    %     if ~isequal(grid.body(body).collocation,collocation{body})
    %         error('oct18_2026_two_drop_image_count_sweep:ChangedGrid', ...
    %             'Collocation nodes changed during the image-count sweep.');
    %     end
    % end
    maxError = drop_image_validate_surface_conditions( ...
        grid,lambda,gamma,U,S1,e,info.options,points,normals);
    residual(k,:) = max(maxError,[],1)./scale;
    [innerProxy(k),outerProxy(k),innerImages(k,:),outerImages(k,:)] = ...
        drop_image_max_source_strengths(grid,lambda,gamma,opt.imageFamilies);
    innerAll(k) = max([innerProxy(k), ...
        innerImages(k,~isnan(innerImages(k,:)))]);
    outerAll(k) = max([outerProxy(k), ...
        outerImages(k,~isnan(outerImages(k,:)))]);
    iterations(k) = info.iterations;
    algebraicResidual(k) = info.originalResidual;
end

%% Report residuals and source-family maxima separately
fprintf('\nIndependent max surface residuals, scaled by |S1|a or |S1|:\n');
disp(table(imageCounts(:),residual(:,1),residual(:,2), ...
    residual(:,3),iterations,algebraicResidual, ...
    'VariableNames',{'NImage','Continuity','NormalVelocity', ...
    'TangentialTraction','GMRESIterations','AlgebraicResidual'}));
familyHeaders = cell(1,nFamilies);
familyLabels = cell(1,nFamilies);
for f = 1:nFamilies
    familyHeaders{f} = char("Image_"+families(f));
    familyLabels{f} = char("image "+families(f));
end
columnNames = [{'NImage','MaxAll','ProxyS'},familyHeaders];
fprintf('Inner source points (exterior-flow coefficients):\n');
disp(array2table([imageCounts(:),innerAll,innerProxy,innerImages], ...
    'VariableNames',columnNames));
fprintf('Outer source points (interior-flow coefficients):\n');
disp(array2table([imageCounts(:),outerAll,outerProxy,outerImages], ...
    'VariableNames',columnNames));
fprintf('All cases use the same collocation and validation nodes.\n');

%% Plot all three physical surface conditions against image count
figure('Name','Two-drop image-count residual sweep', ...
    'Color','w','Position',[100 100 1400 440]);
tiledlayout(1,3,'TileSpacing','compact');
condition = {'velocity continuity','normal velocity', ...
    'tangential traction'};
for j = 1:3
    ax = nexttile;
    semilogy(ax,imageCounts,residual(:,j),'-o','LineWidth',1.4);
    set(ax,'XGrid','on','YGrid','on');
    xlabel(ax,'Image nodes per line');
    ylabel(ax,'Maximum scaled residual');
    title(ax,condition{j});
end
sgtitle('Two drops: fixed geometry, proxy grid, and validation grid');

%% Plot proxy and per-family image maxima for both source locations
figure('Name','Two-drop image-count source magnitudes', ...
    'Color','w','Position',[100 100 1400 520]);
tiledlayout(1,2,'TileSpacing','compact');
strengths = {[innerAll,innerProxy,innerImages], ...
    [outerAll,outerProxy,outerImages]};
location = {'Inner source points (exterior flow)', ...
    'Outer source points (interior flow)'};
legendNames = [{'max all','proxy S'},familyLabels];
for side = 1:2
    ax = nexttile;
    hold(ax,'on');
    for j = 1:2+nFamilies
        semilogy(ax,imageCounts,strengths{side}(:,j),'-o', ...
            'LineWidth',1.2);
    end
    set(ax,'XGrid','on','YGrid','on');
    xlabel(ax,'Image nodes per line');
    ylabel(ax,'Maximum source-strength norm');
    title(ax,location{side});
    legend(ax,legendNames,'Location','best');
end
sgtitle('Proxy and image strengths by source family');
