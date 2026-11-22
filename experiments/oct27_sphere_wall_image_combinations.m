%OCT27_SPHERE_WALL_IMAGE_COMBINATIONS Fixed sphere/wall resistance experiment.
% Deliberate execution: 187 direct SVD solves in split mode, plus two fitting-resolution checks.
% Standard proxy Stokeslets are present in EVERY solve. All wall corrections
% come from StokesHalfspaceFMM; only interior-sphere strengths are fitted.
% T selects T_n, T_n2 and T_t together; T_t and T_nT_n2 can be selected alone.
% Components and coefficients are Cartesian and interleaved by node.

%% FIXED PARAMETERS -- N is exactly 700 and is never swept or rounded.
radius = 1;
mu = 1;
N = 700;
delta_wall = 1e-3*radius;
proxyRadius = 0.63*radius;
wallHeight = 0;
translationVelocity = [1;0;-1];
angularVelocity = [0;0;0];
n_images = 0:5:30;
imageOffsetFraction = 0.05;
svdCutoff = 1e-12;                 % Relative; strict sigma > cutoff*sigma_max.
coarseCollocationRequest = 2*N;    % Existing design: actual count 1353.
capAngle = pi/5;
capPointCount = 1201;
capLayers = 36;
capExponent = 4;
capMinAzimuth = 8;
capAzimuthOversampling = 1;
validationResolution = [83,96];
validationCapAngles = [pi/5,pi/30];
validationCapCounts = [6001,6001];
batchSize = 256;
wallGridCount = 41;
wallGridHalfWidth = 2*radius;
wallRadialCount = 24;
wallAzimuthCount = 48;
wallMinimumRadius = 1e-4*radius;
wallMaximumRadius = radius;
wallCheckImageCounts = [15,30];
wallCheckCombinations = ["SRD","SRTD"];
wallResidualTarget = 1e-10;
runResolutionChecks = true;       % Compare errors on denser grids; refit two cases.
refinedValidationResolution = [127,144];
refinedValidationCapCounts = [12001,12001];
refinedCoarseCollocationRequest = 2*coarseCollocationRequest;
refinedCapPointCount = 2401;
refinedCapLayers = 54;
validationChangeTolerance = 0.05;
fittingChangeTolerance = 0.10;
resolutionResidualFloor = 1e-10;
stressletCombinationMode = "split"; % "split": 31 subsets; "grouped": original 15.
% "custom" uses these strings, e.g. compare the two stresslet groups directly.
customCombinationStrings = ["S","T_t","T_nT_n2","ST_tD","ST_nT_n2D","SRTD"];
outputPrefix = 'oct27_sphere_wall';

%% PATHS AND CONVENTIONS
projectRoot = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(projectRoot,'setpath.m'));
outputDirectory = fullfile(projectRoot,'visuals');
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
if N ~= 700
    error('oct27:FixedN','This experiment requires exactly N=700 proxy sources.');
end
kernelOptions = struct('method','blake','mu',mu,'wall_height',wallHeight);
if stressletCombinationMode == "grouped"
    combinationStrings = ["S","R","T","D","SR","ST","SD","RT", ...
        "RD","TD","SRT","SRD","STD","RTD","SRTD"];
elseif stressletCombinationMode == "split"
    familyTokens = ["S","R","T_nT_n2","T_t","D"];
    combinationStrings = strings(1,31);
    for mask = 1:31
        for family = 1:5
            if bitget(mask,family)
                combinationStrings(mask) = combinationStrings(mask)+familyTokens(family);
            end
        end
        combinationStrings(mask) = replace(combinationStrings(mask),"T_nT_n2T_t","T");
    end
elseif stressletCombinationMode == "custom"
    combinationStrings = customCombinationStrings;
else
    error('oct27:CombinationMode','Use grouped, split or custom combination mode.');
end
selectionMasks = false(numel(combinationStrings),6);
for k = 1:numel(combinationStrings)
    selectionMasks(k,:) = wall_image_combination(combinationStrings(k));
end
if size(unique(selectionMasks,'rows'),1) ~= numel(combinationStrings)
    error('oct27:Combinations','Every requested combination must select a unique basis.');
end
kernelPath = which('blake.stresslet');
if isempty(kernelPath) || isempty(which('blake.potential_dipole'))
    error('oct27:MissingKernels', ...
        'Add the complete StokesHalfspaceFMM/matlab package to the MATLAB path.');
end
kernelFiles = {'blake.velocity','blake.rotlet','blake.stresslet', ...
    'blake.potential_dipole','blake.internal.singularities'};
kernelProvenance = strings(numel(kernelFiles),3);
for k = 1:numel(kernelFiles)
    file = which(kernelFiles{k}); fileInfo = dir(file);
    kernelProvenance(k,:) = [string(kernelFiles{k}),string(file), ...
        string(fileInfo.date)];
end
centre = [0,0,wallHeight+radius+delta_wall];
rigidMotion = [translationVelocity;angularVelocity];
lineTangent = [0,0,-1];            % Along the line, normal to the wall.
lineNormal1 = [0,-1,0];            % cross(lineTangent,e_x), tangent to wall.
lineNormal2 = [-1,0,0];            % cross(lineTangent,lineNormal1).
stressletDirections = [lineNormal1;lineNormal2;lineTangent];
blockLabels = ["S","R","T_n","T_n2","T_t","D"];
blockFamilies = ["S","R","T","T","T","D"];
packageStrengthMultipliers = [1,1,-1,-1,-1,-2];
% R/T strengths have force*length units; D has force*length^2 units.
blockScales = [1,radius,radius,radius,radius,radius^2];
lengthPowers = [0,1,1,1,1,2];
sphere = struct('axes',radius*ones(3,1),'Q',eye(3),'c',centre.');
contactPoint = centre.'+radius*lineTangent.';

%% EXACTLY 700 STANDARD PROXY STOKESLETS AND A FIXED FITTING GRID
% get_sphdesign(700) returns 686, so reuse the existing Fibonacci patch over
% the full proxy sphere instead. This preserves exactly 700 unknown vectors
% without adding a geometry generator or padding/duplicating design nodes.
proxySphere = sphere; proxySphere.axes = proxyRadius*ones(3,1);
proxyOptions = struct('collocationPatchMode',"fibonacci", ...
    'collocationPatchPointCount',N,'collocationPatchAngle',pi);
proxyPoints = collocation_surface_patch(proxySphere, ...
    centre.'+[0;0;proxyRadius],proxyOptions);
if size(proxyPoints,1) ~= N || size(unique(proxyPoints,'rows'),1) ~= N
    error('oct27:ProxyCount','The proxy grid must contain exactly 700 unique points.');
end
[coarseDirections,~] = get_sphdesign(coarseCollocationRequest);
capOptions = struct('collocationPatchMode',"exponential", ...
    'collocationPatchAngle',capAngle,'collocationPatchPointCount',capPointCount, ...
    'collocationPatchLayers',capLayers,'collocationExponent',capExponent, ...
    'collocationMinAzimuth',capMinAzimuth, ...
    'collocationAzimuthOversampling',capAzimuthOversampling);
[capPoints,capMetadata] = collocation_surface_patch(sphere,contactPoint,capOptions);
collocationPoints = uniquetol([centre+radius*coarseDirections;capPoints], ...
    1e-12,'ByRows',true,'DataScale',radius);
collocationRhs = getKmat(collocationPoints,centre)*rigidMotion;
fitVelocityScale = max(vecnorm(reshape(collocationRhs,3,[]).',2,2));
largestColumns = 3*N+18*max(n_images);
if 3*size(collocationPoints,1) < 2*largestColumns
    error('oct27:Oversampling','Fixed fitting grid must have at least two rows per unknown.');
end

%% ONE INDEPENDENT VALIDATION SET, INCLUDING TWO WALL-FACING CAPS
surface = setupsurfquad(ellipsoid_param(radius,radius,radius),validationResolution);
validationPoints = centre+surface.x.';
for k = 1:numel(validationCapAngles)
    options = struct('collocationPatchMode',"fibonacci", ...
        'collocationPatchAngle',validationCapAngles(k), ...
        'collocationPatchPointCount',validationCapCounts(k));
    patch = collocation_surface_patch(sphere,contactPoint,options);
    validationPoints = [validationPoints;patch]; %#ok<AGROW>
end
validationPoints = uniquetol([validationPoints;contactPoint.'], ...
    1e-12,'ByRows',true,'DataScale',radius);
validationDesired = reshape(getKmat(validationPoints,centre)*rigidMotion,3,[]).';
velocityScale = max(vecnorm(validationDesired,2,2));
if velocityScale == 0
    error('oct27:Motion','The prescribed rigid motion must have nonzero speed.');
end

%% INDEPENDENT WALL TARGETS -- NEVER USED IN A FIT
wallCoordinate = linspace(-wallGridHalfWidth,wallGridHalfWidth,wallGridCount);
[wallX,wallY] = meshgrid(wallCoordinate,wallCoordinate);
wallPoints = [wallX(:),wallY(:),wallHeight*ones(numel(wallX),1)];
wallRadii = logspace(log10(wallMinimumRadius),log10(wallMaximumRadius),wallRadialCount);
wallAngles = 2*pi*(0:wallAzimuthCount-1)/wallAzimuthCount;
[wallTheta,wallRho] = meshgrid(wallAngles,wallRadii);
wallPoints = unique([wallPoints;wallRho(:).*cos(wallTheta(:)), ...
    wallRho(:).*sin(wallTheta(:)),wallHeight*ones(numel(wallRho),1)],'rows');

%% GEOMETRY -- THE REFLECTED PAIR GAP IS 2*delta_wall
% Reflection followed by sphere inversion: F(r)=a^2/(2*h-r).
% Its interior fixed point is r*=h-sqrt(h^2-a^2). The equal reflected
% sphere has separation 2*h, so image_line_points sees gap 2*delta_wall.
height = radius+delta_wall;
reflectedCentre = centre-[0,0,2*height];
accumulationRadius = radius^2/(height+sqrt(delta_wall*(2*radius+delta_wall)));
lineInnerRadius = proxyRadius*(1+imageOffsetFraction);
lineEndpoints = [centre+lineInnerRadius*lineTangent; ...
    centre+accumulationRadius*lineTangent];
linePoints = cell(numel(n_images),1);
lineMetadata = cell(numel(n_images),1);
for imageIndex = 1:numel(n_images)
    [linePoints{imageIndex},lineMetadata{imageIndex}] = image_line_points( ...
        centre,reflectedCentre,radius,proxyRadius,n_images(imageIndex), ...
        imageOffsetFraction,radius);
    line = linePoints{imageIndex};
    if abs(lineMetadata{imageIndex}.accumulationRadius-accumulationRadius) > 1e-12*radius ...
            || any(vecnorm(line-centre,2,2) >= radius) ...
            || any(line(:,3) <= wallHeight)
        error('oct27:LineGeometry','Image line failed its sphere/wall geometry check.');
    end
end
fprintf('\nFixed N=%d, delta_wall/a=%g, Rp/a=%g, V=(%g,%g,%g), Omega=(%g,%g,%g).\n', ...
    N,delta_wall/radius,proxyRadius/radius,translationVelocity,angularVelocity);
fprintf('Fitting %d targets; validation %d; largest system %d x %d.\n', ...
    size(collocationPoints,1),size(validationPoints,1), ...
    3*size(collocationPoints,1),largestColumns);
fprintf('Line radii [%g,%g]; accumulation height above wall %.12g.\n', ...
    lineInnerRadius,accumulationRadius,height-accumulationRadius);

%% FLAT RESULTS TABLE AND SAVED COEFFICIENTS
combinationCount = numel(combinationStrings);
caseCount = combinationCount*numel(n_images);
results = table(repmat(combinationStrings(:),numel(n_images),1), ...
    repelem(n_images(:),combinationCount),'VariableNames',{'Combination','ImagePoints'});
metricNames = ["ImageScalars","ScalarUnknowns","Rows","NumericalRank", ...
    "MatrixMiB","SvdBudgetMiB","FittingMaxRelative","FittingRelativeTwoNorm", ...
    "SurfaceMaxRelative","SurfaceMaxAbsolute","CoefficientNorm", ...
    "MaxScalarCoefficient","MaxSourceMagnitude","MaxPackageSourceMagnitude", ...
    "SolveSeconds","AssemblySeconds","ValidationSeconds"];
for name = metricNames, results.(name) = nan(caseCount,1); end
results.BaselineReused = false(caseCount,1);
strengths = cell(caseCount,1);
physicalCoefficients = cell(caseCount,1);
packageCoefficients = cell(caseCount,1);
basisBySolve = cell(caseCount,1);
retainedSingularValues = cell(caseCount,1);
surfacePointResiduals = cell(caseCount,1);
familyMaxima = nan(caseCount,7);    % proxy S, S, R, T_n, T_n2, T_t, D.
wallChecks = table('Size',[0,7], ...
    'VariableTypes',{'string','double','string','double','double','double','logical'}, ...
    'VariableNames',{'Combination','ImagePoints','Block','MaxRelativeWallSpeed', ...
    'RoundoffEstimate','ConstituentScale','MeetsAbsoluteTarget'});
template = struct('family',"S",'points',proxyPoints, ...
    'normalDirections',zeros(0,3),'scale',1);
proxyBasis = template;
timer = tic;
proxyMatrix = wall_image_velocity_matrix(proxyBasis,collocationPoints,kernelOptions);
proxyAssemblySeconds = toc(timer);
factorizationCount = 0;

%% SOURCE COMBINATIONS, WITH ONE SHARED ZERO-IMAGE BASELINE
for imageIndex = 1:numel(n_images)
    imageCount = n_images(imageIndex);
    imageBasis = repmat(template,6,1);
    imageMatrix = zeros(size(proxyMatrix,1),0);
    imageAssemblySeconds = 0;
    if imageCount > 0
        for b = 1:6
            imageBasis(b).family = blockFamilies(b);
            imageBasis(b).points = linePoints{imageIndex};
            imageBasis(b).scale = blockScales(b);
            if blockFamilies(b) == "T"
                imageBasis(b).normalDirections = repmat(stressletDirections(b-2,:),imageCount,1);
            end
        end
        timer = tic;
        imageMatrix = wall_image_velocity_matrix(imageBasis,collocationPoints,kernelOptions);
        imageAssemblySeconds = toc(timer);
    end
    for combinationIndex = 1:combinationCount
        row = (imageIndex-1)*combinationCount+combinationIndex;
        combination = combinationStrings(combinationIndex);
        if imageCount == 0 && combinationIndex > 1
            results(row,:) = results(1,:);
            results.Combination(row) = combination;
            results.BaselineReused(row) = true;
            strengths{row} = strengths{1};
            physicalCoefficients{row} = physicalCoefficients{1};
            packageCoefficients{row} = packageCoefficients{1};
            basisBySolve{row} = basisBySolve{1};
            retainedSingularValues{row} = retainedSingularValues{1};
            surfacePointResiduals{row} = surfacePointResiduals{1};
            familyMaxima(row,:) = familyMaxima(1,:);
            continue
        end
        selected = selectionMasks(combinationIndex,:);
        if imageCount == 0, selected(:) = false; end
        selectedBlocks = find(selected);
        selectedColumns = zeros(1,3*imageCount*numel(selectedBlocks));
        for k = 1:numel(selectedBlocks)
            selectedColumns((k-1)*3*imageCount+(1:3*imageCount)) = ...
                (selectedBlocks(k)-1)*3*imageCount+(1:3*imageCount);
        end
        basis = [proxyBasis;imageBasis(selectedBlocks)];
        matrix = [proxyMatrix,imageMatrix(:,selectedColumns)];
        label = sprintf('%s, n_images=%d',combination,imageCount);
        [strength,d,fitted] = solve_wall_image_system(matrix,collocationRhs,svdCutoff,label);
        factorizationCount = factorizationCount+1;
        clear matrix
        coefficients = image_source_coefficients(basis,strength);
        packageCoeff = coefficients;
        activeMultipliers = [1,packageStrengthMultipliers(selectedBlocks)];
        activePowers = [0,lengthPowers(selectedBlocks)];
        activeIndices = [1,1+selectedBlocks];
        for b = 1:numel(basis)
            familyMaxima(row,activeIndices(b)) = ...
                max(vecnorm(coefficients{b},2,2))/radius^activePowers(b);
            packageCoeff{b} = activeMultipliers(b)*coefficients{b};
        end
        timer = tic;
        pointResidual = zeros(size(validationPoints,1),1);
        for first = 1:batchSize:size(validationPoints,1)
            nodes = first:min(first+batchSize-1,size(validationPoints,1));
            velocity = direct_wall_image_velocity(basis,coefficients, ...
                validationPoints(nodes,:),kernelOptions);
            pointResidual(nodes) = vecnorm(velocity-validationDesired(nodes,:),2,2);
        end
        validationSeconds = toc(timer);
        fittingResidual = reshape(fitted-collocationRhs,3,[]).';
        maxPackageStrength = 0;
        for b = 1:numel(basis)
            maxPackageStrength = max(maxPackageStrength, ...
                max(vecnorm(packageCoeff{b},2,2))/radius^activePowers(b));
        end
        results.ImageScalars(row) = numel(strength)-3*N;
        results.ScalarUnknowns(row) = numel(strength);
        results.Rows(row) = d.rows;
        results.NumericalRank(row) = d.rank;
        results.MatrixMiB(row) = d.matrixMiB;
        results.SvdBudgetMiB(row) = d.workingMiB;
        results.FittingMaxRelative(row) = max(vecnorm(fittingResidual,2,2))/fitVelocityScale;
        results.FittingRelativeTwoNorm(row) = d.fittingRelativeTwoNorm;
        results.SurfaceMaxAbsolute(row) = max(pointResidual);
        results.SurfaceMaxRelative(row) = max(pointResidual)/velocityScale;
        results.CoefficientNorm(row) = d.coefficientNorm;
        results.MaxScalarCoefficient(row) = d.maxScalarCoefficient;
        results.MaxSourceMagnitude(row) = max(familyMaxima(row,:),[],'omitnan');
        results.MaxPackageSourceMagnitude(row) = maxPackageStrength;
        results.SolveSeconds(row) = d.solveSeconds;
        results.AssemblySeconds(row) = proxyAssemblySeconds+imageAssemblySeconds;
        results.ValidationSeconds(row) = validationSeconds;
        strengths{row} = strength;
        physicalCoefficients{row} = coefficients;
        packageCoefficients{row} = packageCoeff;
        basisBySolve{row} = basis;
        retainedSingularValues{row} = d.singularValues;
        surfacePointResiduals{row} = pointResidual/velocityScale;
        fprintf('  Independent surface max %.3e; fit max %.3e; max source %.3e.\n', ...
            results.SurfaceMaxRelative(row),results.FittingMaxRelative(row), ...
            results.MaxSourceMagnitude(row));

        % Independent no slip: total field and every selected source block.
        if ismember(combination,wallCheckCombinations) && ismember(imageCount,wallCheckImageCounts)
            names = ["total","proxy S",blockLabels(selectedBlocks)];
            for check = 1:numel(names)
                wallVelocity = zeros(size(wallPoints));
                wallScale = zeros(size(wallPoints,1),1);
                if check == 1, blocks = 1:numel(basis); else, blocks = check-1; end
                for first = 1:batchSize:size(wallPoints,1)
                    nodes = first:min(first+batchSize-1,size(wallPoints,1));
                    [wallVelocity(nodes,:),wallScale(nodes)] = direct_wall_image_velocity( ...
                        basis(blocks),coefficients(blocks),wallPoints(nodes,:),kernelOptions);
                end
                wallError = max(vecnorm(wallVelocity,2,2))/velocityScale;
                roundoffEstimate = 300*eps*max(wallScale)/velocityScale;
                wallChecks(end+1,:) = {combination,imageCount,names(check),wallError, ...
                    roundoffEstimate,max(wallScale),wallError <= wallResidualTarget}; %#ok<SAGROW>
                if wallError > max(wallResidualTarget,roundoffEstimate)
                    warning('oct27:WallResidual','%s %s wall residual %.3e exceeds roundoff estimate %.3e.', ...
                        label,names(check),wallError,roundoffEstimate);
                end
            end
            fprintf('  Independent wall total residual %.3e.\n', ...
                wallChecks.MaxRelativeWallSpeed(end-numel(basis)));
        end
    end
end
clear imageMatrix proxyMatrix
expectedFactorizations = 1+combinationCount*(numel(n_images)-1);
if factorizationCount ~= expectedFactorizations
    error('oct27:SolveCount','Unexpected number of production factorizations.');
end

%% RESOLUTION CHECKS -- FIX N, GAP, MOTION AND ALL IMAGE GEOMETRY
if runResolutionChecks
    surface = setupsurfquad(ellipsoid_param(radius,radius,radius),refinedValidationResolution);
    refinedValidationPoints = centre+surface.x.';
    for k = 1:numel(validationCapAngles)
        options = struct('collocationPatchMode',"fibonacci", ...
            'collocationPatchAngle',validationCapAngles(k), ...
            'collocationPatchPointCount',refinedValidationCapCounts(k));
        refinedValidationPoints = [refinedValidationPoints; ...
            collocation_surface_patch(sphere,contactPoint,options)]; %#ok<AGROW>
    end
    refinedValidationPoints = uniquetol([refinedValidationPoints;contactPoint.'], ...
        1e-12,'ByRows',true,'DataScale',radius);
    refinedDesired = reshape(getKmat(refinedValidationPoints,centre)*rigidMotion,3,[]).';
    refinedVelocityScale = max(vecnorm(refinedDesired,2,2));
    resolutionCheckRows = unique([1;find(ismember(results.Combination,["SRD","SRTD"]) ...
        & results.ImagePoints==max(n_images))],'stable');
    refinedErrors = zeros(numel(resolutionCheckRows),1);
    for k = 1:numel(resolutionCheckRows)
        row = resolutionCheckRows(k);
        for first = 1:batchSize:size(refinedValidationPoints,1)
            nodes = first:min(first+batchSize-1,size(refinedValidationPoints,1));
            velocity = direct_wall_image_velocity(basisBySolve{row},physicalCoefficients{row}, ...
                refinedValidationPoints(nodes,:),kernelOptions);
            refinedErrors(k) = max(refinedErrors(k), ...
                max(vecnorm(velocity-refinedDesired(nodes,:),2,2))/refinedVelocityScale);
        end
    end
    originalErrors = results.SurfaceMaxRelative(resolutionCheckRows);
    relativeChanges = abs(refinedErrors-originalErrors)./max(originalErrors,resolutionResidualFloor);
    validationResolutionChecks = table(results.Combination(resolutionCheckRows),results.ImagePoints(resolutionCheckRows), ...
        originalErrors,refinedErrors,relativeChanges,relativeChanges<=validationChangeTolerance, ...
        'VariableNames',{'Combination','ImagePoints','OriginalMax','RefinedMax','ScaledChange','Passed'});
    [directions,~] = get_sphdesign(refinedCoarseCollocationRequest);
    options = capOptions;
    options.collocationPatchPointCount = refinedCapPointCount;
    options.collocationPatchLayers = refinedCapLayers;
    refinedCollocationPoints = uniquetol([centre+radius*directions; ...
        collocation_surface_patch(sphere,contactPoint,options)], ...
        1e-12,'ByRows',true,'DataScale',radius);
    refinedRhs = getKmat(refinedCollocationPoints,centre)*rigidMotion;
    resolutionCheckRows = unique([1;find(results.Combination=="SRTD" ...
        & results.ImagePoints==max(n_images))],'stable');
    refinedErrors = zeros(numel(resolutionCheckRows),1);
    refinedSolveSeconds = zeros(numel(resolutionCheckRows),1);
    for k = 1:numel(resolutionCheckRows)
        row = resolutionCheckRows(k);
        basis = basisBySolve{row};
        matrix = wall_image_velocity_matrix(basis,refinedCollocationPoints,kernelOptions);
        [z,d] = solve_wall_image_system(matrix,refinedRhs,svdCutoff, ...
            sprintf('fitting-resolution check %s/%d',results.Combination(row),results.ImagePoints(row)));
        clear matrix
        refinedSolveSeconds(k) = d.solveSeconds;
        coefficients = image_source_coefficients(basis,z);
        for first = 1:batchSize:size(validationPoints,1)
            nodes = first:min(first+batchSize-1,size(validationPoints,1));
            velocity = direct_wall_image_velocity(basis,coefficients, ...
                validationPoints(nodes,:),kernelOptions);
            refinedErrors(k) = max(refinedErrors(k), ...
                max(vecnorm(velocity-validationDesired(nodes,:),2,2))/velocityScale);
        end
    end
    originalErrors = results.SurfaceMaxRelative(resolutionCheckRows);
    relativeChanges = abs(refinedErrors-originalErrors)./max(originalErrors,resolutionResidualFloor);
    fittingResolutionChecks = table(results.Combination(resolutionCheckRows),results.ImagePoints(resolutionCheckRows), ...
        originalErrors,refinedErrors,relativeChanges,relativeChanges<=fittingChangeTolerance, ...
        refinedSolveSeconds,'VariableNames',{'Combination','ImagePoints','OriginalMax', ...
        'RefinedFitMax','ScaledChange','Passed','SolveSeconds'});
    if any(~validationResolutionChecks.Passed) || any(~fittingResolutionChecks.Passed)
        warning('oct27:ResolutionCheck', ...
            'Some resolution checks failed. Inspect saved tables before interpreting the curves.');
    end
else
    validationResolutionChecks = table(); fittingResolutionChecks = table(); %#ok<UNRCH> User can disable resolution checks above.
    refinedValidationPoints = zeros(0,3); refinedCollocationPoints = zeros(0,3);
end

%% TABLES AND NUMERICAL OUTPUT -- SAVE BEFORE FIGURE EXPORT
unknownCounts = array2table(reshape(results.ScalarUnknowns,combinationCount,[]), ...
    'VariableNames',cellstr("n_"+string(n_images)));
unknownCounts = addvars(unknownCounts,combinationStrings(:),'Before',1,'NewVariableNames','Combination');
parameters = struct('radius',radius,'mu',mu,'N',N,'delta_wall',delta_wall, ...
    'proxyRadius',proxyRadius,'wallHeight',wallHeight,'rigidMotion',rigidMotion, ...
    'n_images',n_images,'imageOffsetFraction',imageOffsetFraction, ...
    'svdCutoff',svdCutoff, ...
    'coarseCollocationRequest',coarseCollocationRequest,'capAngle',capAngle, ...
    'capPointCount',capPointCount,'capLayers',capLayers,'capExponent',capExponent, ...
    'capMinAzimuth',capMinAzimuth,'capAzimuthOversampling',capAzimuthOversampling, ...
    'validationResolution',validationResolution,'validationCapAngles',validationCapAngles, ...
    'validationCapCounts',validationCapCounts,'batchSize',batchSize, ...
    'wallGridCount',wallGridCount,'wallGridHalfWidth',wallGridHalfWidth, ...
    'wallRadialCount',wallRadialCount,'wallAzimuthCount',wallAzimuthCount, ...
    'wallMinimumRadius',wallMinimumRadius,'wallMaximumRadius',wallMaximumRadius, ...
    'wallCheckImageCounts',wallCheckImageCounts,'wallCheckCombinations',wallCheckCombinations, ...
    'wallResidualTarget',wallResidualTarget,'runResolutionChecks',runResolutionChecks, ...
    'refinedValidationResolution',refinedValidationResolution, ...
    'refinedValidationCapCounts',refinedValidationCapCounts, ...
    'refinedCoarseCollocationRequest',refinedCoarseCollocationRequest, ...
    'refinedCapPointCount',refinedCapPointCount,'refinedCapLayers',refinedCapLayers, ...
    'validationChangeTolerance',validationChangeTolerance, ...
    'fittingChangeTolerance',fittingChangeTolerance,'resolutionResidualFloor',resolutionResidualFloor, ...
    'actualCollocationCount',size(collocationPoints,1), ...
    'actualValidationCount',size(validationPoints,1),'velocityScale',velocityScale, ...
    'factorizationCount',factorizationCount,'stressletCombinationMode',stressletCombinationMode, ...
    'customCombinationStrings',customCombinationStrings,'matlabVersion',version);
stem = fullfile(outputDirectory,outputPrefix);
writetable(results,[stem '_image_combinations.csv']);
writetable(unknownCounts,[stem '_unknown_counts.csv']);
writetable(wallChecks,[stem '_wall_checks.csv']);
if runResolutionChecks
    writetable(validationResolutionChecks,[stem '_validation_resolution_checks.csv']);
    writetable(fittingResolutionChecks,[stem '_fitting_resolution_checks.csv']);
end
save([stem '_image_combinations.mat'],'parameters','results','unknownCounts', ...
    'combinationStrings','selectionMasks','blockLabels','stressletDirections','packageStrengthMultipliers', ...
    'blockScales','lengthPowers','centre','proxyPoints','collocationPoints', ...
    'capMetadata','validationPoints','wallPoints','linePoints','lineMetadata', ...
    'lineEndpoints','accumulationRadius','strengths','physicalCoefficients', ...
    'packageCoefficients','basisBySolve','retainedSingularValues','surfacePointResiduals', ...
    'familyMaxima','wallChecks','validationResolutionChecks','fittingResolutionChecks', ...
    'refinedValidationPoints','refinedCollocationPoints','kernelProvenance','-v7.3');

%% PLOT SAVED RESULTS AGAINST IMAGE UNKNOWNS, THEN SHOW THE NODES
% This function loads only saved metrics and coordinates. It never solves.
% Rerun this section, or call the plotting function with the MAT filename,
% to regenerate both figures without repeating any experiment calculations.
[figureHandle,nodeFigureHandle] = plot_oct27_sphere_wall_image_combinations( ...
    [stem '_image_combinations.mat'],stem);
residualCurves = reshape(results.SurfaceMaxRelative,combinationCount,[]);
bestResidual = min(residualCurves,[],2);
disp(table(combinationStrings(:),bestResidual,'VariableNames',{'Combination','BestSurfaceResidual'}));
disp(unknownCounts);
if runResolutionChecks, disp(validationResolutionChecks); disp(fittingResolutionChecks); end
fprintf('Saved %d cases from %d production factorizations, with N=%d fixed, under %s.\n', ...
    caseCount,factorizationCount,N,outputDirectory);
