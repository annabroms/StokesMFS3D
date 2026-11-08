%% Exterior-only stagnant-cap bubble: compare cone source families
% Purpose: compare cone Stokeslets with Stokeslets plus three oriented
% stresslet blocks at one nondegenerate cap angle. The uniform interior
% proxy, cone sites, force, and SVD tolerance are identical in both cases.
% Only exterior velocity and traction are solved; no interior fluid enters.
%
% Edit: phi, n_level, ringLevels, the source and collocation azimuth counts,
% collocAngularExtent, and the clustering controls below. The default
% has five cone levels and 10 off-edge collocation rings PER region. Edit
% familySets, derivativeDirections, and use_orth_dir_images to compare any
% supported cone blocks and orientation convention.
% The script raises collocation azimuths as families add unknowns; it reports
% the actual count and a dense SVD memory preflight before each solve.
% Toggle useColumnScaling for unit-norm columns, and useAreaWeights for
% square-root local surface-area weighting of the boundary equations.
%
% Checks: use the same shifted validation grids for every case. Report
% rank, force balance, transverse drift, torque, speed error against
% Sadhal-Johnson (53), and four independently normalized boundary maxima.
% Repeat at higher density and closer to the edge, flagging unstable maxima.
%
% Output: experiments/results/oct24_stagnant_cap_bubble_residuals_families.csv
% and case-labelled geometry, residual, and source-strength PNGs. Existing
% files with those names are overwritten when this script is run.
% Run: run('drops/experiments/oct24_stagnant_cap_bubble_residuals.m')
clear; close all;
experimentDirectory = fileparts(mfilename('fullpath'));
outputDirectory = fullfile(experimentDirectory,'results');
if ~isfolder(outputDirectory), mkdir(outputDirectory); end

% CLUSTERING CONTROLS (sources and collocation are independent):
% sourceClustering / collocClustering: "geometric" or "tapered".
% For N levels ordered from outermost to closest, k=0:N-1:
%   geometric: offsets = extent*ratio.^k (smaller ratio reaches closer).
%   tapered: offsets = extent*exp(-sigma*(sqrt(N)-sqrt(N-k))).
% Larger sigma reaches closer; sigma=2 to 4 is a starting range to investigate.
% Only *ClusterRatio is used in geometric mode; only *ClusterSigma in tapered.
% Source extent = innerSourceStartDistance (length); N = n_level.
% Collocation extent = collocAngularExtent (radians); N = ringLevels per side.
% Keep ringLevels >= 2*n_level. Increasing N in tapered mode reaches closer
% AND fills intermediate scales; adding source families does not add scales.
% Example source switch: sourceClustering="tapered", sourceClusterSigma=4.
% Collocation need not use the same mode or sigma. For a first comparison,
% keep it fixed; geometric collocClusterRatio=0.5 fills scales more densely
% than 0.1, which puts ten rings from 0.2 down to only 2e-10 radians.
% Alternatively choose collocClustering="tapered" and tune collocClusterSigma.
% Check printed minima and offsets, not just the ratio: N differs on each side.
% Retain angular coverage of the residual band; do not only add closer rings.
% Very shallow sources require finer azimuths in BOTH n_az_source/n_az_colloc;
% column scaling and surface-area weighting do not supply missing resolution.
% SOURCE-DEPTH REFINEMENT (Stokeslets only in familySets below):
% useSourceRefinement adds local collocation patches at the closest source
% level(s). For source depth d, radialFactors*d is the meridional distance
% from the edge, on BOTH sides; azimuthFactors*d is the arc offset around
% each source's surface projection. This resolves narrow source footprints
% without uniformly increasing every ring to O(a/d) azimuths. The global
% grids and all source sites remain fixed. sourceRefinementLevels selects
% how many closest source levels are covered; factors control patch spacing.
% Shifted, denser copies are added to independent validation as well.
commonInput = struct('radius',1,'muOut',1,'phi',pi/4, ...
    'innerProxyRatio',0.63,'sourceRequest',700, ...
    'collocationRatio',2,'svdTolerance',1e-8, ...
    'useColumnScaling',false,'useAreaWeights',false, ...
    'gravity',[0 0 -1],'densityDifference',1, ...
    'n_level',5,'ringLevels',10, ...
    'innerSourceStartDistance',0.35, ...
    'useSourceRefinement',true,'sourceRefinementLevels',1, ...
    'sourceRefinementRadialFactors',[0.25 0.5 1 2], ...
    'sourceRefinementAzimuthFactors',[-2 -1 0 1 2], ...
    'sourceClustering',"tapered",'sourceClusterSigma',3, ...
    'sourceClusterRatio',0.1,'n_az_source',64, ...
    'n_az_colloc',80,'collocAngularExtent',0.5, ...
    'collocClustering',"geometric",'collocClusterSigma',6, ...
    'collocClusterRatio',0.5,'use_orth_dir_images',false);
familySets = {"S"}; % Use Stokeslets only at refining source nodes for now.
derivativeDirections = "none";
assert(numel(familySets)==numel(derivativeDirections));
caseCount = numel(familySets);
caseOptions = cell(caseCount,1);
geometries = cell(caseCount,1);
caseNames = strings(caseCount,1);
caseTags = strings(caseCount,1);
for k = 1:caseCount
    input = commonInput;
    input.refinementFamilies = familySets{k};
    input.derivative_direction = derivativeDirections(k);
    caseOptions{k} = bubble_options(input);
    geometries{k} = bubble_geometry(caseOptions{k});
    caseNames(k) = strjoin(caseOptions{k}.refinementFamilies,'+');
    if derivativeDirections(k)~="none"
        caseNames(k) = caseNames(k)+"+S_y("+derivativeDirections(k)+")";
    end
    caseTags(k) = regexprep(caseNames(k),'[^A-Za-z0-9]+','_');
    if k>1
        assert(isequal(geometries{k}.sourcePoints, ...
            geometries{1}.sourcePoints) && ...
            isequal(geometries{k}.ringOffsets, ...
            geometries{1}.ringOffsets), ...
            'Family cases must share the cone sites and ring offsets.');
        if geometries{k}.coneFamilyCount> ...
                geometries{1}.coneFamilyCount
            assert(geometries{k}.collocationCount> ...
                geometries{1}.collocationCount, ...
                'Additional cone blocks must increase collocation points.');
        end
    end
end
assert(numel(unique(caseTags))==caseCount, ...
    'Source-family cases need distinct output tags.');
referenceOptions = caseOptions{1};
Fg = (4*pi*referenceOptions.radius^3/3)* ...
    referenceOptions.densityDifference*referenceOptions.gravity(:);
theoreticalUz = Fg(3)/stagnant_cap_exact_drag( ...
    referenceOptions.phi,0,referenceOptions.muOut, ...
    referenceOptions.radius);
maxCollocation = max(cellfun(@(g) g.collocationCount,geometries));
finestOffset = geometries{1}.ringOffsets(end);
validation = bubble_validation_grid( ...
    referenceOptions,maxCollocation,finestOffset);
fineOptions = referenceOptions;
fineOptions.validationBaseRequest = max( ...
    fineOptions.validationBaseRequest,2*maxCollocation);
fineOptions = bubble_options(fineOptions);
fineOffset = min(finestOffset/2, ...
    validation.minimumEdgeDistance/referenceOptions.radius);
fineValidation = bubble_validation_grid( ...
    fineOptions,2*maxCollocation,fineOffset);
fprintf(['Shared validation: %d points at %.5g a from edge; ', ...
    'finer %d points at %.5g a.\n'], ...
    validation.count,validation.minimumEdgeDistance/referenceOptions.radius, ...
    fineValidation.count, ...
    fineValidation.minimumEdgeDistance/referenceOptions.radius);

sourceSites = zeros(caseCount,1);
useColumnScaling = false(caseCount,1);
useAreaWeights = false(caseCount,1);
coneSites = zeros(caseCount,1);
coneFamilyCount = zeros(caseCount,1);
coefficientCount = zeros(caseCount,1);
ringLevels = zeros(caseCount,1);
ringAzimuthCount = zeros(caseCount,1);
collocationCount = zeros(caseCount,1);
svdRows = zeros(caseCount,1);
svdColumns = zeros(caseCount,1);
retainedRank = zeros(caseCount,1);
knownDenseGB = zeros(caseCount,1);
minimumSourceDistanceOverA = zeros(caseCount,1);
minimumAddedRingDistanceOverA = zeros(caseCount,1);
minimumCollocationEdgeDistanceOverA = zeros(caseCount,1);
numericalUz = zeros(caseCount,1);
relativeSpeedError = zeros(caseCount,1);
forceRelativeResidual = zeros(caseCount,1);
transverseDrift = zeros(caseCount,1);
normalizedTorque = zeros(caseCount,1);
validationUnresolved = false(caseCount,1);
boundaryMaxima = zeros(caseCount,4);
validationChange = zeros(caseCount,4);
peakDistanceOverA = zeros(caseCount,4);
for k = 1:caseCount
    opt = caseOptions{k};
    geometry = geometries{k};
    fprintf('\nCase %s: %d cone families, %d rings per region, ', ...
        caseNames(k),geometry.coneFamilyCount,opt.ringLevels);
    fprintf('%d collocation azimuths.\n',geometry.ringAzimuthCount);
    fprintf(['Closest source %.6g a from interface; added ring ', ...
        '%.6g a and actual collocation %.6g a from edge.\n'], ...
        geometry.minimumSourceDistance/opt.radius, ...
        geometry.minimumAddedEdgeDistance/opt.radius, ...
        geometry.minimumEdgeDistance/opt.radius);
    fprintf('Source clustering %s (ratio %.4g, sigma %.4g).\n', ...
        opt.sourceClustering,opt.sourceClusterRatio,opt.sourceClusterSigma);
    if opt.n_level>0
        sourceOffsets = stagnant_cap_cluster_offsets(opt.n_level, ...
            opt.innerSourceStartDistance,opt.sourceClustering, ...
            opt.sourceClusterRatio,opt.sourceClusterSigma);
        fprintf('Source depths/a:');
        fprintf(' %.5g',sourceOffsets/opt.radius);
        fprintf('\n');
    end
    fprintf('Collocation clustering %s (ratio %.4g, sigma %.4g); offsets/rad:', ...
        opt.collocClustering,opt.collocClusterRatio,opt.collocClusterSigma);
    fprintf(' %.5g',geometry.ringOffsets);
    fprintf('\n');
    fprintf('Source-depth refinement: %d added local collocation points.\n', ...
        geometry.sourceRefinementCount);
    solution = bubble_solve(opt,Fg,geometry);
    check = bubble_check_validation( ...
        solution,opt,validation,fineValidation);
    sourceSites(k) = geometry.sourceCount;
    useColumnScaling(k) = opt.useColumnScaling;
    useAreaWeights(k) = opt.useAreaWeights;
    coneSites(k) = geometry.coneSourceCount;
    coneFamilyCount(k) = geometry.coneFamilyCount;
    coefficientCount(k) = geometry.coefficientCount;
    ringLevels(k) = geometry.collocationRingCount(1);
    ringAzimuthCount(k) = geometry.ringAzimuthCount;
    collocationCount(k) = geometry.collocationCount;
    svdRows(k) = solution.matrixSize(1);
    svdColumns(k) = solution.matrixSize(2);
    retainedRank(k) = solution.retainedRank;
    knownDenseGB(k) = solution.svdKnownDenseMB/1000;
    minimumSourceDistanceOverA(k) = ...
        geometry.minimumSourceDistance/opt.radius;
    minimumAddedRingDistanceOverA(k) = ...
        geometry.minimumAddedEdgeDistance/opt.radius;
    minimumCollocationEdgeDistanceOverA(k) = ...
        geometry.minimumEdgeDistance/opt.radius;
    numericalUz(k) = solution.U(3);
    relativeSpeedError(k) = abs(solution.U(3)-theoreticalUz)/ ...
        abs(theoreticalUz);
    forceRelativeResidual(k) = solution.forceResidual;
    transverseDrift(k) = solution.transverseDrift;
    normalizedTorque(k) = norm(solution.torque)/ ...
        (norm(Fg)*opt.radius);
    validationUnresolved(k) = check.unresolved;
    boundaryMaxima(k,:) = check.report.maxima;
    validationChange(k,:) = check.relativeChange;
    fprintf(['Uz %.8g, theory %.8g, relative error %.3g, ', ...
        'rank %d/%d, force residual %.3g, torque %.3g.\n'], ...
        numericalUz(k),theoreticalUz,relativeSpeedError(k), ...
        retainedRank(k),min(solution.matrixSize), ...
        forceRelativeResidual(k),normalizedTorque(k));
    [geometry3D,geometryMeridional] = bubble_plot_geometry( ...
        geometry,opt,'on');
    residualFigure = bubble_plot_residuals( ...
        check.report,validation,solution,opt,'on');
    [collocationFigure,collocationReport] = bubble_plot_collocation_residuals( ...
        solution,opt,'on');
    [distanceFigure,residualProfiles] = bubble_plot_residuals_distance( ...
        check.report,validation,solution,opt,'on');
    peakDistanceOverA(k,:) = [residualProfiles.peakDistanceOverA];
    for condition = 1:4
        fprintf('%s: peak at s/a %.6g (%.3g source depths), azimuth %.6g rad.\n', ...
            check.report.conditionNames(condition), ...
            residualProfiles(condition).peakDistanceOverA, ...
            residualProfiles(condition).peakDistanceOverA* ...
                opt.radius/geometry.minimumSourceDistance, ...
            residualProfiles(condition).peakAzimuth);
    end
    % Retain pointwise fields so plots can be regenerated without a new solve.
    save(fullfile(outputDirectory,sprintf( ...
        'oct24_stagnant_cap_bubble_residuals_%s_data.mat',char(caseTags(k)))), ...
        'solution','opt','check','residualProfiles','collocationReport');
    [strength3D,strengthMeridional] = bubble_plot_strengths( ...
        solution,opt,'on');
    handles = [geometry3D,geometryMeridional,residualFigure, ...
        collocationFigure,distanceFigure,strength3D,strengthMeridional];
    suffixes = {'geometry_3d','geometry_meridional','residuals', ...
        'residuals_collocation','residuals_distance', ...
        'strengths_3d','strengths_meridional'};
    for j = 1:numel(handles)
        exportgraphics(handles(j),fullfile(outputDirectory, ...
            sprintf('oct24_stagnant_cap_bubble_residuals_%s_%s.png', ...
            char(caseTags(k)),suffixes{j})),'Resolution',180);
    end
end

theoreticalSpeed = repmat(theoreticalUz,caseCount,1);
validationCount = repmat(validation.count,caseCount,1);
fineValidationCount = repmat(fineValidation.count,caseCount,1);
minimumValidationEdgeDistanceOverA = repmat( ...
    validation.minimumEdgeDistance/referenceOptions.radius,caseCount,1);
minimumFineValidationEdgeDistanceOverA = repmat( ...
    fineValidation.minimumEdgeDistance/referenceOptions.radius,caseCount,1);
summary = table(caseNames,useColumnScaling,useAreaWeights, ...
    sourceSites,coneSites,coneFamilyCount, ...
    coefficientCount,ringLevels,ringAzimuthCount,collocationCount, ...
    svdRows,svdColumns,retainedRank,knownDenseGB, ...
    minimumSourceDistanceOverA,minimumAddedRingDistanceOverA, ...
    minimumCollocationEdgeDistanceOverA, ...
    validationCount,fineValidationCount, ...
    minimumValidationEdgeDistanceOverA, ...
    minimumFineValidationEdgeDistanceOverA, ...
    numericalUz,theoreticalSpeed, ...
    relativeSpeedError,forceRelativeResidual,transverseDrift, ...
    normalizedTorque,validationUnresolved);
summary.capNormalMax = boundaryMaxima(:,1);
summary.mobileNormalMax = boundaryMaxima(:,2);
summary.capTangentMax = boundaryMaxima(:,3);
summary.mobileTractionMax = boundaryMaxima(:,4);
summary.capNormalValidationChange = validationChange(:,1);
summary.mobileNormalValidationChange = validationChange(:,2);
summary.capTangentValidationChange = validationChange(:,3);
summary.mobileTractionValidationChange = validationChange(:,4);
summary.capNormalPeakDistanceOverA = peakDistanceOverA(:,1);
summary.mobileNormalPeakDistanceOverA = peakDistanceOverA(:,2);
summary.capTangentPeakDistanceOverA = peakDistanceOverA(:,3);
summary.mobileTractionPeakDistanceOverA = peakDistanceOverA(:,4);
% Save shared placement parameters so runs can be compared without ambiguity.
for name = ["sourceClustering","sourceClusterRatio","sourceClusterSigma", ...
        "collocClustering","collocClusterRatio","collocClusterSigma", ...
        "innerSourceStartDistance","collocAngularExtent","n_az_source", ...
        "useSourceRefinement","sourceRefinementLevels"]
    summary.(name) = repmat(referenceOptions.(name),caseCount,1);
end
summary.sourceRefinementRadialFactors = repmat( ...
    string(mat2str(referenceOptions.sourceRefinementRadialFactors)),caseCount,1);
summary.sourceRefinementAzimuthFactors = repmat( ...
    string(mat2str(referenceOptions.sourceRefinementAzimuthFactors)),caseCount,1);
writetable(summary,fullfile(outputDirectory, ...
    'oct24_stagnant_cap_bubble_residuals_families.csv'));
disp(summary);
