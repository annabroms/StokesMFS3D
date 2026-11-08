%% Compare exterior bubble formulations on identical Stokeslet geometry
% Both methods use force projection and the same known-force completion.
% Only cap velocity rows differ: split normal/tangents versus full Cartesian.
% Run explicitly after setpath. Both modes and matched resolution checks use
% one configured SVD cutoff. Edit formulationInput
% before running to override any settings below. No source coefficients are
% compared. All fields use common nonzero physical reference scales.
% Example small configuration (still an experiment, not an automatic test):
% formulationInput.options = struct('sourceRequest',50,'n_level',1, ...
%     'ringLevels',2,'n_az_source',6,'n_az_colloc',8, ...
%     'useSourceRefinement',false);
% Revised outputs use the oct29_bubble_cap_formulations_force_completion prefix.
% Existing weighted-force experiment outputs remain available.
if ~exist('formulationInput','var'), formulationInput = struct(); end
if isfield(formulationInput,'forceWeight')
    fprintf('Removing obsolete forceWeight: both cap modes now use algebraic force completion.\n');
    formulationInput = rmfield(formulationInput,'forceWeight');
end
if isfield(formulationInput,'runSensitivity')
    formulationInput = rmfield(formulationInput,'runSensitivity');
end
defaults = struct('options',struct(), ...
    'runResolutionCheck',true,'resolutionFactor',1.25, ...
    'visibility','on','saveOutputs',true,'planeCount',161,'profileAzimuthCount',4, ...
    'outputDirectory',fullfile(fileparts(mfilename('fullpath')),'results'));
assert(isstruct(formulationInput) && isscalar(formulationInput));
assert(all(ismember(fieldnames(formulationInput),fieldnames(defaults))));
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(formulationInput,names{k})
        formulationInput.(names{k}) = defaults.(names{k});
    end
end
assert(islogical(formulationInput.runResolutionCheck) && isscalar(formulationInput.runResolutionCheck));
assert(islogical(formulationInput.saveOutputs) && isscalar(formulationInput.saveOutputs));
assert(ismember(string(formulationInput.visibility),["on","off"]));
validateattributes(formulationInput.resolutionFactor,{'numeric'},{'scalar','finite','>',1});
validateattributes(formulationInput.planeCount,{'numeric'},{'scalar','integer','>=',21});
validateattributes(formulationInput.profileAzimuthCount,{'numeric'},{'scalar','integer','positive'});
commonInput = struct('radius',1,'muOut',1,'phi',pi/4, ...
    'innerProxyRatio',0.63,'sourceRequest',700,'collocationRatio',2, ...
    'svdTolerance',1e-8,'maxSvdMemoryGB',8, ...
    'useColumnScaling',false,'useAreaWeights',false, ...
    'conditionWeights',[1 1 1],'useConditionCountScaling',false, ...
    'gravity',[0 0 -1],'densityDifference',1,'n_level',5,'ringLevels',10, ...
    'innerSourceStartDistance',0.35,'sourceClustering',"tapered", ...
    'sourceClusterSigma',3,'sourceClusterRatio',0.1, ...
    'n_az_source',64,'n_az_colloc',80,'collocAngularExtent',0.5, ...
    'collocClustering',"geometric",'collocClusterRatio',0.5, ...
    'useSourceRefinement',true,'sourceRefinementLevels',1, ...
    'sourceRefinementRadialFactors',[0.25 0.5 1 2], ...
    'sourceRefinementAzimuthFactors',[-2 -1 0 1 2], ...
    'refinementFamilies',"S",'derivative_direction',"none");
names = fieldnames(formulationInput.options);
for k = 1:numel(names)
    commonInput.(names{k}) = formulationInput.options.(names{k});
end
if commonInput.n_level==0, commonInput.useSourceRefinement = false; end
opt = bubble_options(commonInput);
assert(opt.phi>0 && opt.phi<pi,'This edge-profile experiment requires 0<phi<pi.');
assert(isequal(opt.refinementFamilies,"S") && opt.derivative_direction=="none");
assert(~opt.useColumnScaling && ~opt.useAreaWeights && ...
    ~opt.useConditionCountScaling && isequal(opt.conditionWeights,[1 1 1]), ...
    'Comparison requires the documented isotropic cap weighting.');
capModes = ["split","cartesian"];
settings = struct('forceTreatment',"projection_completion",'capModes',capModes, ...
    'completionScale',opt.completionScale);
geometry = bubble_geometry(opt);
caseNames = "baseline";
caseOptions = {opt};
geometries = {geometry};
if formulationInput.runResolutionCheck
    next = opt;
    next.sourceRequest = ceil(formulationInput.resolutionFactor*opt.sourceRequest);
    next.ringLevels = max(opt.ringLevels+1,ceil(formulationInput.resolutionFactor*opt.ringLevels));
    if opt.n_level>0
        next.n_level = opt.n_level+1;
        next.ringLevels = max(next.ringLevels,2*next.n_level);
        next.n_az_source = ceil(sqrt(formulationInput.resolutionFactor)*opt.n_az_source);
        next.n_az_colloc = ceil(sqrt(formulationInput.resolutionFactor)*opt.n_az_colloc);
    end
    next = bubble_options(next);
    caseOptions{end+1} = next;
    geometries{end+1} = bubble_geometry(next);
    caseNames(end+1) = "refined";
end
caseCount = numel(caseOptions);
methodNames = ["split_cap","cartesian_cap"];
fprintf('All experiment solves use relative SVD cutoff %.6g.\n',opt.svdTolerance);
% Check every configured geometry before spending time on the first dense solve.
for c = 1:caseCount
    g = geometries{c};
    assert(caseOptions{c}.svdTolerance==opt.svdTolerance, ...
        'Every experiment case must use the same SVD cutoff.');
    for method = 1:2
        estimate = bubble_svd_memory(3*g.collocationCount,g.coefficientCount);
        fprintf('%s/%s planned matrix %d x %d; estimated dense SVD %.4g GB.\n', ...
            caseNames(c),methodNames(method),estimate.rows,estimate.columns,estimate.estimatedGB);
        assert(estimate.estimatedGB<=caseOptions{c}.maxSvdMemoryGB, ...
            'Configured %s/%s exceeds the %.4g GB budget; reduce resolution or raise the budget.', ...
            caseNames(c),methodNames(method),caseOptions{c}.maxSvdMemoryGB);
    end
end
Fg = (4*pi*opt.radius^3/3)*opt.densityDifference*opt.gravity(:);
zeta = stagnant_cap_exact_drag(opt.phi,0,opt.muOut,opt.radius);
theoryU = Fg/zeta;
Ustar = norm(theoryU);
tractionStar = opt.muOut*Ustar/opt.radius;
forceStar = norm(Fg);
assert(Ustar>0 && tractionStar>0 && forceStar>0);
maxCollocation = max(cellfun(@(g) g.collocationCount,geometries));
allCollocation = zeros(0,3);
finestOffset = opt.edgeBandAngle/2;
for k = 1:caseCount
    allCollocation = [allCollocation;geometries{k}.points]; %#ok<AGROW>
    if ~isempty(geometries{k}.ringOffsets)
        finestOffset = min(finestOffset,min(geometries{k}.ringOffsets));
    end
end

%% Shared independent surface grids: rotated base, shifted rings, no edge
% Six times the largest M leaves room to remove accidental coincidences.
auditOptions = opt;
auditOptions.validationMultiplier = max(6,opt.validationMultiplier);
auditOptions = bubble_options(auditOptions);
validation = bubble_validation_grid(auditOptions,maxCollocation,finestOffset);
fineOffset = min(finestOffset/2,validation.minimumEdgeDistance/opt.radius);
fineValidation = bubble_validation_grid(auditOptions,2*maxCollocation,fineOffset);
grids = {validation,fineValidation};
gridNames = ["primary","fine"];
for k = 1:2
    v = grids{k};
    coincidences = ismembertol(v.points,allCollocation,1e-11*opt.radius, ...
        'ByRows',true,'DataScale',1);
    keep = ~coincidences & v.edgeDistance>1e-12*opt.radius;
    v.baseCount = nnz(keep(1:v.baseCount));
    for name = ["points","normals","theta","cap","mobile","edgeDistance","nearEdge"]
        v.(name) = v.(name)(keep,:);
    end
    v.count = size(v.points,1);
    v.minimumEdgeDistance = min(v.edgeDistance);
    assert(v.count>=5*k*maxCollocation && all(v.edgeDistance>0));
    assert(~any(ismembertol(v.points,allCollocation,1e-11*opt.radius, ...
        'ByRows',true,'DataScale',1)));
    assert(any(v.cap) && any(v.mobile));
    fprintf(['%s validation: %d points, excluded %d coincidences/edge points; ', ...
        'minimum arc distance %.8g (%.8g a); cap/mobile %.8g/%.8g a.\n'], ...
        gridNames(k),v.count,nnz(~keep),v.minimumEdgeDistance, ...
        v.minimumEdgeDistance/opt.radius,min(v.edgeDistance(v.cap))/opt.radius, ...
        min(v.edgeDistance(v.mobile))/opt.radius);
    grids{k} = v;
end
validation = grids{1};
fineValidation = grids{2};
assert(fineValidation.count>validation.count && ...
    fineValidation.minimumEdgeDistance<validation.minimumEdgeDistance);

%% Matched solves and independently measured errors
solutions = cell(caseCount,2);
reports = cell(caseCount,2,2);
summaryRows = cell(2*caseCount,1);
for c = 1:caseCount
    current = caseOptions{c};
    g = geometries{c};
    fprintf('\n%s: unchanged source law %s; collocation law %s.\n', ...
        caseNames(c),current.sourceClustering,current.collocClustering);
    fprintf(['Physical scales U*=%.8g, t*=%.8g, F*=%.8g; ', ...
        'shared completion beta %.8g; force balance imposed algebraically.\n'], ...
        Ustar,tractionStar,forceStar,current.completionScale);
    for method = 1:2
        s = bubble_solve(current,Fg,g,[],capModes(method));
        assert(isequal(s.grid.sourcePoints,g.sourcePoints) && ...
            isequal(s.grid.points,g.points) && isequal(s.grid.normals,g.normals) && ...
            isequal(s.grid.cap,g.cap) && isequal(s.grid.mobile,g.mobile));
        assert(abs(s.Uref/Ustar-1)<1e-13 && abs(s.tractionScale/tractionStar-1)<1e-13);
        solutions{c,method} = s;
        for gridIndex = 1:2
            reports{c,method,gridIndex} = bubble_validate(s,grids{gridIndex},current,true);
        end
        primary = reports{c,method,1};
        fine = reports{c,method,2};
        first = [primary.maxCapVelocity,primary.maxima([2 4])];
        second = [fine.maxCapVelocity,fine.maxima([2 4])];
        change = abs(second-first)./max(max(abs(first),abs(second)),1e-12);
        row = struct('Case',caseNames(c),'Method',methodNames(method), ...
            'Radius',opt.radius,'MuOut',opt.muOut,'Phi',opt.phi, ...
            'FgX',Fg(1),'FgY',Fg(2),'FgZ',Fg(3), ...
            'Ux',s.U(1),'Uy',s.U(2),'Uz',s.U(3),'TheoryUz',theoryU(3), ...
            'RelativeTranslationError',norm(s.U-theoryU)/Ustar, ...
            'TransverseDrift',s.transverseDrift, ...
            'CapVelocity',second(1),'MobileNormal',second(2),'MobileTraction',second(3), ...
            'CapNormal',fine.maxima(1),'CapTangent',fine.maxima(3), ...
            'PrimaryCapVelocity',first(1),'PrimaryMobileNormal',first(2), ...
            'PrimaryMobileTraction',first(3),'ForceResidual',s.forceResidual, ...
            'ForceDefectX',s.FH(1)+Fg(1),'ForceDefectY',s.FH(2)+Fg(2), ...
            'ForceDefectZ',s.FH(3)+Fg(3),'NormalizedTorque',norm(s.torque)/(forceStar*opt.radius), ...
            'BaseSources',g.baseSourceCount,'ConeSources',g.coneSourceCount, ...
            'Sources',g.sourceCount,'StrengthUnknowns',g.coefficientCount, ...
            'TotalUnknowns',s.matrixSize(2),'Collocation',g.collocationCount, ...
            'CapTargets',g.capCount,'MobileTargets',g.mobileCount, ...
            'BaseTargets',g.baseCollocationCount,'PatchTargets',g.sourceRefinementCount, ...
            'RingTargets',nnz(g.collocationRingIndex>0), ...
            'MatrixRows',s.matrixSize(1),'MatrixColumns',s.matrixSize(2), ...
            'Rank',s.retainedRank,'RelativeCutoff',s.svdTolerance, ...
            'AbsoluteCutoff',s.svdAbsoluteThreshold,'SigmaMax',s.sigmaMax, ...
            'KnownDenseGB',s.svdKnownDenseMB/1000,'EstimatedGB',s.svdMemoryEstimateGB, ...
            'Ustar',Ustar,'TractionStar',tractionStar,'ForceStar',forceStar, ...
            'CompletionScale',current.completionScale,'CapMode',s.capMode, ...
            'ForceFreeResidual',s.forceFreeResidual, ...
            'CompletionForceResidual',s.completionForceResidual, ...
            'AreaWeights',s.useAreaWeights,'ColumnScaling',s.useColumnScaling, ...
            'ConditionCountScaling',s.useConditionCountScaling, ...
            'NormalWeight',s.conditionWeights(1),'CapTangentWeight',s.conditionWeights(2), ...
            'MobileTractionWeight',s.conditionWeights(3), ...
            'PrimaryValidation',validation.count,'FineValidation',fineValidation.count, ...
            'MinimumEdgeOverA',fineValidation.minimumEdgeDistance/opt.radius, ...
            'MinimumCollocationEdgeOverA',g.minimumEdgeDistance/opt.radius, ...
            'MinimumSourceDepthOverA',g.minimumSourceDistance/opt.radius, ...
            'CapValidationChange',change(1),'NormalValidationChange',change(2), ...
            'TractionValidationChange',change(3),'ValidationUnresolved',any(change>0.1));
        summaryRows{2*(c-1)+method} = row;
        fprintf(['%s/%s: U=[% .8g % .8g % .8g]; theory error %.4g; ', ...
            'fine cap/normal/traction %.4g/%.4g/%.4g; force %.4g.\n'], ...
            caseNames(c),methodNames(method),s.U,row.RelativeTranslationError, ...
            second,s.forceResidual);
    end
    assert(isequal(solutions{c,1}.basis,solutions{c,2}.basis));
end
summaryTable = struct2table(vertcat(summaryRows{:}));
disp(summaryTable(:,{'Case','Method','Uz','RelativeTranslationError', ...
    'CapVelocity','MobileNormal','MobileTraction','ForceResidual','Rank'}));

%% Physical field differences; coefficient agreement is not meaningful
comparisonRows = cell(2*caseCount,1);
for c = 1:caseCount
    for gridIndex = 1:2
        v = grids{gridIndex};
        old = reports{c,1,gridIndex};
        new = reports{c,2,gridIndex};
        du = old.fields.velocity-new.fields.velocity;
        dr = old.fields.relativeVelocity-new.fields.relativeVelocity;
        dt = old.fields.tangentTraction-new.fields.tangentTraction;
        duNorm = vecnorm(du,2,2)/Ustar;
        dtNorm = vecnorm(dt,2,2)/tractionStar;
        relativeNorm = vecnorm(dr,2,2)/Ustar;
        normalDifference = abs(sum(dr.*v.normals,2))/Ustar;
        deltaU = norm(solutions{c,1}.U-solutions{c,2}.U)/Ustar;
        speedErrorSum = norm(solutions{c,1}.U-theoryU)/Ustar+ ...
            norm(solutions{c,2}.U-theoryU)/Ustar;
        speedRatio = 0;
        if speedErrorSum>0, speedRatio = deltaU/speedErrorSum; end
        assert(deltaU<=speedErrorSum+1e-12);
        capBound = old.maxCapVelocity+new.maxCapVelocity;
        normalBound = old.maxima(2)+new.maxima(2);
        tractionBound = old.maxima(4)+new.maxima(4);
        assert(max(relativeNorm(v.cap))<=capBound+1e-10*max(1,capBound));
        assert(max(normalDifference(v.mobile))<=normalBound+1e-10*max(1,normalBound));
        assert(max(dtNorm(v.mobile))<=tractionBound+1e-10*max(1,tractionBound));
        comparisonRows{2*(c-1)+gridIndex} = struct('Case',caseNames(c), ...
            'Grid',gridNames(gridIndex),'Count',v.count, ...
            'MinimumEdgeOverA',v.minimumEdgeDistance/opt.radius, ...
            'DeltaTranslation',deltaU,'TranslationErrorSum',speedErrorSum, ...
            'TranslationConsistencyRatio',speedRatio, ...
            'DeltaVelocity',max(duNorm),'DeltaRelativeVelocity',max(relativeNorm), ...
            'DeltaTangentialTraction',max(dtNorm), ...
            'DeltaCapTraction',max(dtNorm(v.cap)), ...
            'DeltaMobileTraction',max(dtNorm(v.mobile)), ...
            'DeltaCapRelativeVelocity',max(relativeNorm(v.cap)), ...
            'CapResidualSum',capBound,'DeltaMobileNormal',max(normalDifference(v.mobile)), ...
            'MobileNormalResidualSum',normalBound,'MobileTractionResidualSum',tractionBound);
    end
end
comparisonTable = struct2table(vertcat(comparisonRows{:}));
disp(comparisonTable(:,{'Case','Grid','DeltaTranslation','DeltaVelocity', ...
    'DeltaTangentialTraction','DeltaCapTraction','DeltaMobileTraction'}));

%% Empirical discretization assessment at one fixed SVD cutoff
% Resolution changes estimate discretization sensitivity, not cutoff sensitivity.
E = nan(2,3);
refined = find(caseNames=="refined");
if ~isempty(refined)
    for method = 1:2
        base = reports{1,method,2}.fields;
        probe = reports{refined,method,2}.fields;
        E(method,:) = [norm(solutions{refined,method}.U-solutions{1,method}.U)/Ustar, ...
            max(vecnorm(probe.velocity-base.velocity,2,2))/Ustar, ...
            max(vecnorm(probe.tangentTraction-base.tangentTraction,2,2))/tractionStar];
    end
end
baselineDifference = [comparisonTable.DeltaTranslation(2), ...
    comparisonTable.DeltaVelocity(2),comparisonTable.DeltaTangentialTraction(2)];
budget = sum(E,1);
agreementRatios = baselineDifference./budget;
agreementRatios(budget==0 & baselineDifference==0) = 0;
agreementRatios(budget==0 & baselineDifference>0) = Inf;
agreementStatus = "provisional: matched resolution check not run";
if formulationInput.runResolutionCheck
    refinedDifference = [comparisonTable.DeltaTranslation(2*refined), ...
        comparisonTable.DeltaVelocity(2*refined), ...
        comparisonTable.DeltaTangentialTraction(2*refined)];
    if any(summaryTable.ValidationUnresolved)
        agreementStatus = "unresolved: independent maxima change by more than 10 percent";
    elseif any(agreementRatios>1)
        agreementStatus = "unresolved: discrepancy exceeds observed resolution sensitivity at fixed cutoff";
    elseif any(refinedDifference>max(max(1.1*baselineDifference,budget),1e-12))
        agreementStatus = "unresolved: refined discrepancy exceeds observed sensitivity";
    else
        agreementStatus = "empirically consistent with measured resolution sensitivity at fixed cutoff";
    end
end
agreementTable = table(["translation";"velocity";"tangential traction"], ...
    baselineDifference.',E(1,:).',E(2,:).', ...
    budget.',agreementRatios.',repmat(agreementStatus,3,1), ...
    'VariableNames',{'Quantity','Difference','SplitRefinement','CartesianRefinement', ...
    'EstimatedBudget','Ratio','Status'});
disp(agreementTable);
fprintf('Agreement: %s. Estimators are empirical; cutoff sensitivity is not assessed.\n',agreementStatus);

%% Surface maps on the fine independent grid, with paired color limits
v = fineValidation;
baseReports = {reports{1,1,2},reports{1,2,2}};
psi = atan2(v.normals(:,2),v.normals(:,1));
eTheta = [cos(v.theta).*cos(psi),cos(v.theta).*sin(psi),-sin(v.theta)];
surfaceValues = cell(2,2);
for method = 1:2
    relative = baseReports{method}.fields.relativeVelocity;
    surfaceValues{1,method} = sum(relative.*eTheta,2)/Ustar;
    surfaceValues{2,method} = vecnorm(relative,2,2)/Ustar;
end
figureHandles = gobjects(0);
figureSuffixes = strings(0);
surfaceLimits = zeros(2,2);
f = figure('Name','Bubble-relative surface velocity','Color','w', ...
    'Visible',formulationInput.visibility,'Position',[40 40 1350 1000]);
layout = tiledlayout(f,2,2,'TileSpacing','compact');
edgePsi = linspace(0,2*pi,240).';
edge = opt.centre+opt.radius*[sin(opt.phi)*cos(edgePsi), ...
    sin(opt.phi)*sin(edgePsi),cos(opt.phi)*ones(size(edgePsi))];
% Raise the edge line slightly for visibility only; field points are unchanged.
displayEdge = opt.centre+(1+1e-3)*(edge-opt.centre);
for quantity = 1:2
    unionValues = [surfaceValues{quantity,1};surfaceValues{quantity,2}];
    upper = max(max(abs(unionValues)),eps);
    if quantity==1
        surfaceLimits(quantity,:) = [-upper upper];
        quantityLabel = '(u-U) dot e_\theta / U_*';
    else
        surfaceLimits(quantity,:) = [0 upper];
        quantityLabel = 'total relative speed ||u-U|| / U_*';
    end
    for method = 1:2
        ax = nexttile(layout,2*(quantity-1)+method);
        scatter3(ax,v.points(:,1),v.points(:,2),v.points(:,3),8, ...
            surfaceValues{quantity,method},'filled');
        hold(ax,'on');
        plot3(ax,displayEdge(:,1),displayEdge(:,2),displayEdge(:,3),'k-','LineWidth',2);
        text(ax,opt.centre(1),opt.centre(2),opt.centre(3)+1.05*opt.radius, ...
            'stagnant cap','HorizontalAlignment','center');
        axis(ax,'equal'); view(ax,35,25); grid(ax,'on');
        clim(ax,surfaceLimits(quantity,:)); colormap(ax,turbo);
        ylabel(colorbar(ax),quantityLabel);
        title(ax,sprintf('%s: %s\nU_z/U_* = %.5g; force residual %.3g', ...
            methodNames(method),quantityLabel,solutions{1,method}.U(3)/Ustar, ...
            solutions{1,method}.forceResidual));
        xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
    end
end
sgtitle(layout,{sprintf('Each recovered translation subtracted; minimum edge distance %.4g a', ...
    v.minimumEdgeDistance/opt.radius), ...
    'Cap: zero relative velocity. Mobile region: tangential slip is allowed.', ...
    'Top: signed polar component (positive towards increasing theta). Bottom: total speed includes normal error.'});
figureHandles(end+1) = f;
figureSuffixes(end+1) = "surface_velocity";

f = figure('Name','Surface velocity differences','Color','w', ...
    'Visible',formulationInput.visibility,'Position',[40 40 1350 600]);
layout = tiledlayout(f,1,2);
surfaceDifference = baseReports{1}.fields.velocity-baseReports{2}.fields.velocity;
relativeDifference = baseReports{1}.fields.relativeVelocity-baseReports{2}.fields.relativeVelocity;
differenceValues = [vecnorm(surfaceDifference,2,2),vecnorm(relativeDifference,2,2)]/Ustar;
differenceLimits = [0,max(max(differenceValues,[],'all'),eps)];
labels = ["Laboratory velocity difference / U_*","Relative velocity difference / U_*"];
for k = 1:2
    ax = nexttile(layout,k);
    scatter3(ax,v.points(:,1),v.points(:,2),v.points(:,3),8,differenceValues(:,k),'filled');
    hold(ax,'on'); plot3(ax,edge(:,1),edge(:,2),edge(:,3),'k-');
    axis(ax,'equal'); view(ax,35,25); grid(ax,'on');
    clim(ax,differenceLimits); colormap(ax,turbo); colorbar(ax);
    title(ax,labels(k)); xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
end
figureHandles(end+1) = f;
figureSuffixes(end+1) = "surface_difference";

%% All imposed boundary conditions, plus the legacy split-cap diagnostics
conditionIndices = [5 2 4 1 3];
conditionMasks = {v.cap,v.mobile,v.mobile,v.cap,v.cap};
conditionValues = cell(5,2);
for method = 1:2
    r = baseReports{method};
    conditionValues(:,method) = {r.point.capVelocity/Ustar; ...
        r.point.normal/Ustar;r.point.mobileTraction/tractionStar; ...
        r.point.normal/Ustar;r.point.capTangent/Ustar};
end
residualLimits = zeros(5,2);
f = figure('Name','Paired boundary residuals','Color','w', ...
    'Visible',formulationInput.visibility,'Position',[40 40 1400 1500]);
layout = tiledlayout(f,5,2,'TileSpacing','compact');
residualAxes = gobjects(5,2);
for condition = 1:5
    mask = conditionMasks{condition};
    values = [conditionValues{condition,1}(mask);conditionValues{condition,2}(mask)];
    colors = log10(max(values,1e-16));
    limits = [min(colors),max(colors)];
    if limits(1)==limits(2), limits = limits+[-0.5 0.5]; end
    residualLimits(condition,:) = limits;
    for method = 1:2
        residualAxes(condition,method) = nexttile(layout,2*(condition-1)+method);
    end
end
for method = 1:2
    plotSettings = struct('axes',residualAxes(:,method), ...
        'conditionIndices',conditionIndices,'colorLimits',residualLimits,'logFloor',1e-16);
    bubble_plot_residuals(baseReports{method},v,solutions{1,method},opt, ...
        formulationInput.visibility,char(methodNames(method)),plotSettings);
    for condition = 1:5
        ax = residualAxes(condition,method);
        assert(isequal(clim(ax),residualLimits(condition,:)));
        title(ax,methodNames(method)+": "+string(ax.Title.String));
    end
end
sgtitle(layout,'Fine independent grid: applicable regions only; common limits across each pair');
figureHandles(end+1) = f;
figureSuffixes(end+1) = "residuals";

%% Common polar profiles and close edge views, with azimuthal diagnostics
minimumAngle = v.minimumEdgeDistance/opt.radius;
offsets = logspace(log10(minimumAngle),log10(opt.edgeBandAngle),100).';
profileTheta = sort(unique([((0:399).'+0.37)*pi/400;opt.phi-offsets;opt.phi+offsets]));
profileTheta = profileTheta(profileTheta>0 & profileTheta<pi & ...
    abs(profileTheta-opt.phi)>1e-12);
profilePsi = 0.217+2*pi*(0:formulationInput.profileAzimuthCount-1)/ ...
    formulationInput.profileAzimuthCount;
[thetaGrid,psiGrid] = ndgrid(profileTheta,profilePsi);
profileNormals = [sin(thetaGrid(:)).*cos(psiGrid(:)), ...
    sin(thetaGrid(:)).*sin(psiGrid(:)),cos(thetaGrid(:))];
profilePoints = opt.centre+opt.radius*profileNormals;
assert(~any(ismembertol(profilePoints,allCollocation,1e-11*opt.radius, ...
    'ByRows',true,'DataScale',1)));
profileTangents = [cos(thetaGrid(:)).*cos(psiGrid(:)), ...
    cos(thetaGrid(:)).*sin(psiGrid(:)),-sin(thetaGrid(:))];
profileVelocity = zeros(numel(profileTheta),numel(profilePsi),2);
for method = 1:2
    field = bubble_evaluate(solutions{1,method},profilePoints,opt);
    profileVelocity(:,:,method) = reshape(sum(field.relativeVelocity.*profileTangents,2)/ ...
        Ustar,size(thetaGrid));
end
profileDifference = profileVelocity(:,:,1)-profileVelocity(:,:,2);
profileMean = squeeze(mean(profileVelocity,2));
profileMeanDifference = mean(profileDifference,2);
% Insert a NaN at the excluded edge so plotted traces do not bridge it.
[plotTheta,permutation] = sort([profileTheta;opt.phi]);
plotVelocity = [profileMean;nan(1,2)];
plotVelocity = plotVelocity(permutation,:);
plotDifference = [profileMeanDifference;NaN];
plotDifference = plotDifference(permutation);
f = figure('Name','Surface tangential velocity profiles','Color','w', ...
    'Visible',formulationInput.visibility,'Position',[40 40 1450 900]);
layout = tiledlayout(f,2,2);
for panel = 1:4
    ax = nexttile(layout,panel); hold(ax,'on');
    if mod(panel,2)==1
        plot(ax,plotTheta,plotVelocity(:,1),'b-','LineWidth',1.4,'DisplayName','split cap mean');
        plot(ax,plotTheta,plotVelocity(:,2),'r--','LineWidth',1.4,'DisplayName','Cartesian cap mean');
        for azimuth = 1:numel(profilePsi)
            for method = 1:2
                trace = [profileVelocity(:,azimuth,method);NaN];
                if method==1, color = [0.65 0.75 1]; else, color = [1 0.7 0.7]; end
                plot(ax,plotTheta,trace(permutation),'Color',color,'HandleVisibility','off');
            end
        end
        ylabel(ax,'(u-U) dot e_\theta / U_*');
        legend(ax,'Location','best');
    else
        plot(ax,plotTheta,plotDifference,'k-','LineWidth',1.4);
        for azimuth = 1:numel(profilePsi)
            trace = [profileDifference(:,azimuth);NaN];
            plot(ax,plotTheta,trace(permutation),'Color',[0.65 0.65 0.65]);
        end
        ylabel(ax,'split minus Cartesian / U_*');
    end
    xline(ax,opt.phi,'k:'); xlabel(ax,'polar angle theta'); grid(ax,'on');
    if panel>2, xlim(ax,opt.phi+[-1 1]*opt.edgeBandAngle); else, xlim(ax,[0 pi]); end
end
sgtitle(layout,sprintf('Shared off-edge profiles at %d azimuths; closest distance %.4g a', ...
    numel(profilePsi),min(abs(profileTheta-opt.phi))));
figureHandles(end+1) = f;
figureSuffixes(end+1) = "profiles";

%% Shared exterior meridional grid: physical laboratory velocity and difference
planeAxis = unique([linspace(-3,3,formulationInput.planeCount), ...
    linspace(-1.15,1.15,formulationInput.planeCount)]);
[planeX,planeZ] = meshgrid(opt.centre(1)+opt.radius*planeAxis, ...
    opt.centre(3)+opt.radius*planeAxis);
exterior = (planeX-opt.centre(1)).^2+(planeZ-opt.centre(3)).^2>=opt.radius^2;
planePoints = [planeX(exterior),opt.centre(2)*ones(nnz(exterior),1),planeZ(exterior)];
planeVelocity = nan([size(planeX),3,2]);
for method = 1:2
    field = bubble_evaluate(solutions{1,method},planePoints,opt);
    for component = 1:3
        values = nan(size(planeX));
        values(exterior) = field.velocity(:,component)/Ustar;
        planeVelocity(:,:,component,method) = values;
    end
end
assert(all(isnan(planeVelocity(repmat(~exterior,1,1,3,2)))));
planeDifference = planeVelocity(:,:,:,1)-planeVelocity(:,:,:,2);
planeSpeed = {sqrt(sum(planeVelocity(:,:,:,1).^2,3)), ...
    sqrt(sum(planeVelocity(:,:,:,2).^2,3)),sqrt(sum(planeDifference.^2,3))};
pairedMaximum = max([planeSpeed{1}(exterior);planeSpeed{2}(exterior)]);
planeLimits = [0,max(pairedMaximum,eps);0,max(pairedMaximum,eps); ...
    0,max(max(planeSpeed{3}(exterior)),eps)];
f = figure('Name','Exterior laboratory velocity comparison','Color','w', ...
    'Visible',formulationInput.visibility,'Position',[40 40 1650 600]);
layout = tiledlayout(f,1,3);
stride = max(1,floor(numel(planeAxis)/25));
sample = 1:stride:numel(planeAxis);
angle = linspace(0,2*pi,240);
for panel = 1:3
    ax = nexttile(layout,panel);
    surf(ax,planeX,planeZ,zeros(size(planeX)),planeSpeed{panel}, ...
        'EdgeColor','none','FaceColor','interp');
    view(ax,2); hold(ax,'on');
    if panel<=2
        vectors = planeVelocity(:,:,:,panel);
        label = methodNames(panel)+": norm(u) / U_*";
    else
        vectors = planeDifference;
        label = "Laboratory difference / U_*";
    end
    arrowScale = 0.18*opt.radius;
    quiver(ax,planeX(sample,sample),planeZ(sample,sample), ...
        arrowScale*vectors(sample,sample,1),arrowScale*vectors(sample,sample,3), ...
        0,'k','LineWidth',0.6);
    fill(ax,opt.centre(1)+opt.radius*cos(angle), ...
        opt.centre(3)+opt.radius*sin(angle),'w','EdgeColor','k');
    clim(ax,planeLimits(panel,:)); colormap(ax,turbo); colorbar(ax);
    axis(ax,'equal'); axis(ax,'tight'); xlabel(ax,'x'); ylabel(ax,'z'); title(ax,label);
end
sgtitle(layout,'Same exterior grid and arrow scale; bubble interior masked');
figureHandles(end+1) = f;
figureSuffixes(end+1) = "meridional";

%% One geometry object for both baseline solves
[f3,f2] = bubble_plot_geometry(geometry,opt,formulationInput.visibility);
geometryCoordinateDifference = max(abs(solutions{1,1}.grid.sourcePoints- ...
    solutions{1,2}.grid.sourcePoints),[],'all');
collocationCoordinateDifference = max(abs(solutions{1,1}.grid.points- ...
    solutions{1,2}.grid.points),[],'all');
assert(geometryCoordinateDifference==0 && collocationCoordinateDifference==0);
sgtitle(f3,sprintf('Split and Cartesian cap modes: identical sources/targets; coordinate differences %.1g / %.1g', ...
    geometryCoordinateDifference,collocationCoordinateDifference));
sgtitle(f2,'Split and Cartesian cap modes share every base/cone source and cap/mobile target');
figureHandles = [figureHandles,f3,f2];
figureSuffixes = [figureSuffixes,"geometry_3d","geometry_meridional"];
plotMetadata = struct('surfaceLimits',surfaceLimits,'differenceLimits',differenceLimits, ...
    'residualLimits',residualLimits,'residualFloor',1e-16,'planeLimits',planeLimits, ...
    'planeAxis',planeAxis,'arrowScale',arrowScale,'profilePsi',profilePsi, ...
    'minimumProfileDistanceOverA',min(abs(profileTheta-opt.phi)), ...
    'sourceCoordinateDifference',geometryCoordinateDifference, ...
    'collocationCoordinateDifference',collocationCoordinateDifference);

%% Save reproducible tables, physical fields, figures and measured conclusions
if formulationInput.saveOutputs
    outputDirectory = formulationInput.outputDirectory;
    if ~isfolder(outputDirectory), mkdir(outputDirectory); end
    prefix = fullfile(outputDirectory,'oct29_bubble_cap_formulations_force_completion');
    writetable(summaryTable,[prefix,'_summary.csv']);
    writetable(comparisonTable,[prefix,'_comparison.csv']);
    writetable(agreementTable,[prefix,'_agreement.csv']);
    validationTables = cell(2,1);
    components = ["X","Y","Z"];
    for gridIndex = 1:2
        v = grids{gridIndex};
        t = table(repmat(gridNames(gridIndex),v.count,1),v.points(:,1),v.points(:,2), ...
            v.points(:,3),v.cap,v.edgeDistance/opt.radius, ...
            'VariableNames',{'Grid','X','Y','Z','Cap','EdgeDistanceOverA'});
        for method = 1:2
            r = reports{1,method,gridIndex};
            for component = 1:3
                t.(methodNames(method)+"Velocity"+components(component)) = r.fields.velocity(:,component);
                t.(methodNames(method)+"TangentialTraction"+components(component)) = ...
                    r.fields.tangentTraction(:,component);
            end
            capResidual = r.point.capVelocity;
            capResidual(~v.cap) = NaN;
            capTangentResidual = r.point.capTangent;
            capTangentResidual(~v.cap) = NaN;
            mobileResidual = r.point.mobileTraction;
            mobileResidual(~v.mobile) = NaN;
            t.(methodNames(method)+"CapVelocityResidual") = capResidual;
            t.(methodNames(method)+"NormalVelocityResidual") = r.point.normal;
            t.(methodNames(method)+"CapTangentResidual") = capTangentResidual;
            t.(methodNames(method)+"TangentialTractionResidual") = mobileResidual;
        end
        validationTables{gridIndex} = t;
    end
    validationTable = vertcat(validationTables{:});
    writetable(validationTable,[prefix,'_validation.csv']);
    profileTable = table(thetaGrid(:),psiGrid(:), ...
        reshape(profileVelocity(:,:,1),[],1),reshape(profileVelocity(:,:,2),[],1), ...
        profileDifference(:),'VariableNames',{'Theta','Azimuth', ...
        'SplitRelativeUthetaOverUstar','CartesianRelativeUthetaOverUstar','DifferenceOverUstar'});
    writetable(profileTable,[prefix,'_profiles.csv']);
    save([prefix,'_data.mat'],'formulationInput','settings','caseOptions','caseNames', ...
        'solutions','geometries','grids','reports','Fg','theoryU','Ustar','tractionStar', ...
        'forceStar','summaryTable','comparisonTable','agreementTable','agreementStatus', ...
        'profileTheta','profilePsi','profileVelocity','profileDifference', ...
        'planeX','planeZ','exterior','planeVelocity','planeDifference','plotMetadata','-v7.3');
    for k = 1:numel(figureHandles)
        exportgraphics(figureHandles(k),prefix+"_"+figureSuffixes(k)+".png",'Resolution',180);
    end
    fid = fopen([prefix,'_results.md'],'w');
    assert(fid>=0,'Cannot open results note.');
    fileCleanup = onCleanup(@() fclose(fid)); % Close the results note on failure too.
    fprintf(fid,'# Exterior bubble formulation comparison\n\n');
    fprintf(fid,'%s. These resolution-sensitivity estimates are empirical, not rigorous bounds.\n\n',agreementStatus);
    fprintf(fid,['Laboratory exterior flow at rest at infinity. Cap velocity is relative ', ...
        'to each recovered translation. All field differences use common references.\n\n']);
    fprintf(fid,'U*=%.8g, traction*=%.8g, force*=%.8g; theoretical Uz=%.8g.\n\n', ...
        Ustar,tractionStar,forceStar,theoryU(3));
    fprintf(fid,['Velocity/traction rows are divided by these references. ', ...
        'Both methods use the same force projection and known-force completion; ', ...
        'there are no fitted force rows. Unit isotropic cap weights; area, ', ...
        'count and column scaling off. Shared completion scale %.6g.\n\n'],opt.completionScale);
    fprintf(fid,['Independent surface grids: %d / %d points; minimum edge arc distances ', ...
        '%.8g a / %.8g a. Exact edge excluded.\n\n'],validation.count,fineValidation.count, ...
        validation.minimumEdgeDistance/opt.radius,fineValidation.minimumEdgeDistance/opt.radius);
    fprintf(fid,'| method | Uz | theory error | cap velocity | mobile normal | mobile traction | force | matrix | rank | absolute cutoff |\n');
    fprintf(fid,'| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: |\n');
    for method = 1:2
        r = summaryTable(method,:);
        fprintf(fid,'| %s | %.8g | %.4g | %.4g | %.4g | %.4g | %.4g | %d x %d | %d | %.4g |\n', ...
            methodNames(method),r.Uz,r.RelativeTranslationError,r.CapVelocity,r.MobileNormal, ...
            r.MobileTraction,r.ForceResidual,r.MatrixRows,r.MatrixColumns,r.Rank,r.AbsoluteCutoff);
    end
    fprintf(fid,'\nTranslation / maximum laboratory velocity / maximum tangential-traction differences: %.6g / %.6g / %.6g.\n\n', ...
        baselineDifference);
    fprintf(fid,'All solves use relative SVD cutoff %.6g; cutoff sensitivity is not assessed.\n\n',opt.svdTolerance);
    fprintf(fid,'Difference / observed resolution-refinement budget: %.6g / %.6g / %.6g.\n\n',agreementRatios);
    fprintf(fid,['See the summary, comparison and agreement CSVs for all cases and raw ', ...
        'force defects, force-free source residuals and completion force residuals. ', ...
        'Cap traction is an unconstrained reaction. Source coefficients are not compared.\n\n']);
    fprintf(fid,['Surface figure: the top row is the signed polar component of u-U ', ...
        '(positive towards increasing theta); the bottom row is total relative speed, ', ...
        'which includes normal-velocity error. Both must vanish on the cap; mobile ', ...
        'tangential slip is allowed. Read speed error and force balance alongside these maps.\n']);
    clear fileCleanup
    fprintf('Saved comparison outputs with prefix %s.\n',prefix);
end
