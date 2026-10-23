%DEMO_TWO_DROPS_WITH_IMAGES Compare image discretisations and optional dense solve.
% Run explicitly from the repository root or after adding repository paths.
clear;close all;
runDenseSVD = false;            % Set true for the full-system reference.
plotBody = 1;                   % Show this drop in the surface-residual plots.
modes = 2;                      % Use 1:2 to compare both image configurations.
% Fibonacci uses two overlapping caps centred on the contact ray. Exponential
% uses one cap with capAngles(1), sum(capPointCounts), and graded polar rings.
capMode = "exponential";          % "fibonacci" or "exponential".
capAngles = [pi/5 pi/60];        % Wide/narrow angular radii (36 and 3 degrees).
capPointCounts = [800 200];      % Wide/narrow counts; 1,000 points in total.
capExponent = 3;                % Exponential only: larger clusters nearer contact.
capLayers = 18;                 % Exponential only: number of polar rings.
% imageFamilies tokens: S Stokeslet, R rotlet, D potential dipole, T legacy
% stresslet, T_t line-parallel, T_n first transverse, T_n2 second transverse
% stresslet, Q_s trace-free axial Stokeslet quadrupole. T excludes T_*.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(root));
q = [-1.005 0 0;1.005 0 0];
a = ones(2,1); S1 = ones(2,1); e = repmat([1 0 0],2,1);
base = struct('N',700,'collocationRatio',2,'nearThreshold',0.15, ...
    'N_image',40,'imageFamilies',"SDRT_nT_t",'proxyRadius',0.63, ... %[] for adaptive choice of N_image.
    'use_orth_dir_images',false, ...
    'outerProxyRatio',2,'outerImageOffsetFraction',0.001, ...
    'interiorViscosity',10,'svdCutoff',1e-12, ...
    'gmresTolerance',1e-11,'maxIterations',350,'batchSize',128, ...
    'validationCount',1000,'validationGapCount',401, ...
    'capMode',capMode,'capAngles',capAngles,'capPointCounts',capPointCounts, ...
    'capExponent',[capExponent capExponent],'capLayers',[capLayers capLayers], ...
    'useRowScaling',true,'plot_grid',1);
modeName = {'interior images','interior and exterior images'};
useImages = [true true];
mirror = [false true];
points = []; normals = [];
plotValue = cell(2,2);
plotPoints = cell(2,2);
scale = max(abs(S1(plotBody)),realmin)*[a(plotBody),a(plotBody),1];
nodeSet = {'solve grid','independent nodes'};
ray = q(3-plotBody,:)-q(plotBody,:);
ray = ray/norm(ray);
for mode = modes
    opt = base; opt.useImages = useImages(mode); opt.mirror_ext = mirror(mode);
    opt.solverMode = "oneBody";
    [Ui,li,gi,g,ii] = solve_drops_with_images(q,a,S1,e,opt);
    if isempty(points)
        [points,normals] = drop_image_validation_grid( ...
            g,opt.validationCount,opt.validationGapCount,max(base.capAngles));
    end
    [ei,pointErrorIter] = drop_image_validate_surface_conditions( ...
        g,li,gi,Ui,S1,e,ii.options,points,normals);
    % Evaluate the same physical conditions at the actual solve nodes.
    solvePoints = {g.body.collocation};
    solveNormals = {g.body.normals};
    [~,pointErrorSolve] = drop_image_validate_surface_conditions( ...
        g,li,gi,Ui,S1,e,ii.options,solvePoints,solveNormals);
    if runDenseSVD
        opt.solverMode = "denseSVD";
        [Ud,ld,gd,~,id] = solve_drops_with_images(q,a,S1,e,opt);
        [ed,~] = drop_image_validate_surface_conditions( ...
            g,ld,gd,Ud,S1,e,id.options,points,normals);
    end
    fprintf('\n%s: A is %d x %d, GMRES iterations %d\n', ...
        modeName{mode},ii.systemSize,ii.iterations);
    for body = 1:2
        fprintf('  body %d: proxy %d, inner images %d, outer images %d, collocation %d\n', ...
            body,g.body(body).N,g.body(body).N_image, ...
            size(g.body(body).exteriorImagePoints,1), ...
            size(g.body(body).collocation,1));
    end
    fprintf('  iterative original scaled residual %.3e\n', ...
        ii.originalResidual);
    fprintf('  independent maxima [continuity normal traction]:\n');
    fprintf('    iterative %.3e %.3e %.3e\n',max(ei,[],1));
    fprintf('  iterative translations:\n');
    disp(Ui);
    if runDenseSVD
        fprintf('  relative translation difference %.3e\n', ...
            norm(Ui-Ud,'fro')/max(norm(Ud,'fro'),eps));
        fprintf('  dense original scaled residual %.3e\n', ...
            id.originalResidual);
        fprintf('    dense     %.3e %.3e %.3e\n',max(ed,[],1));
        fprintf('  dense translations:\n');
        disp(Ud);
    end
    % Retain only the selected body's nodes and scaled log10 residuals.
    plotPoints{mode,1} = solvePoints{plotBody};
    plotPoints{mode,2} = points{plotBody};
    plotValue{mode,1} = log10(max(pointErrorSolve{plotBody}./scale,eps));
    plotValue{mode,2} = log10(max(pointErrorIter{plotBody}./scale,eps));
    % Locate continuity peaks relative to the contact ray for cap tuning.
    for row = 1:2
        [value,node] = max(plotValue{mode,row}(:,1));
        n = (plotPoints{mode,row}(node,:)-q(plotBody,:))/a(plotBody);
        theta = acosd(max(-1,min(1,n*ray.')));
        fprintf('  drop %d %s: max log10 continuity %.2f at %.2f degrees from contact\n', ...
            plotBody,nodeSet{row},value,theta);
    end
    drop_image_plot_source_strengths( ...
        g,li,gi,sprintf('%s (iterative)',modeName{mode}));
end

% Match each condition's colour range across both grids and all solved modes.
values = vertcat(plotValue{modes,:});
limits = [min(values,[],1);max(values,[],1)];
constant = limits(1,:)==limits(2,:);
limits(2,constant) = limits(2,constant)+1;
condition = {'continuity','normal velocity','tangential traction'};
for mode = modes
    figure('Name','Two-drop image MFS surface residuals','Color','w', ...
        'Position',[100 100 1400 800]);
    tiledlayout(2,3,'TileSpacing','compact');
    for row = 1:2
        x = plotPoints{mode,row};
        for k = 1:3
            ax = nexttile;
            scatter3(ax,x(:,1),x(:,2),x(:,3),8,plotValue{mode,row}(:,k),'filled');
            axis(ax,'equal'); view(ax,3);
            set(ax,'XGrid','on','YGrid','on','ZGrid','on');
            clim(ax,limits(:,k).');
            title(ax,sprintf('%s: %s',nodeSet{row},condition{k}));
            cb = colorbar(ax);
            cb.Label.String = 'log_{10} scaled residual';
        end
    end
    sgtitle(sprintf('Drop %d: %s (iterative)',plotBody,modeName{mode}));
    colormap(turbo);
end
