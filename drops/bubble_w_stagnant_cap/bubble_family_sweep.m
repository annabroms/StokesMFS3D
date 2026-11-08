function result = bubble_family_sweep(input,controls)
%BUBBLE_FAMILY_SWEEP Compare enhancing bases on identical fit/audit points.
% Preserve all source sites, force completion, scaling and SVD cutoff.
% More families mean more columns; this is a fixed-site, not fixed-DOF test.
if nargin<2, controls = struct(); end
defaults = struct('familySets',{{"S","SD","SR","ST_t","ST_n", ...
    "ST_n2","SQ_s","ST_tT_nT_n2","ST_tT_nT_n2R", ...
    "SDRT_tT_nT_n2Q_s"}},'derivativeDirections',strings(0,1), ...
    'maxMemoryGB',8,'maxSvdRows',30000,'maxSvdColumns',10000, ...
    'maxCollocationPoints',10000,'planOnly',false,'makePlots',true, ...
    'figureVisible','on','outputDirectory','');
names = fieldnames(controls);
assert(all(ismember(names,fieldnames(defaults))),'Unknown sweep control.');
for k = 1:numel(names), defaults.(names{k}) = controls.(names{k}); end
c = defaults;
assert(iscell(c.familySets) && ~isempty(c.familySets));
for name = ["maxSvdRows","maxSvdColumns","maxCollocationPoints"]
    validateattributes(c.(name),{'numeric'},{'scalar','integer','positive','finite'});
end
validateattributes(c.maxMemoryGB,{'numeric'},{'scalar','positive','finite'});
assert(islogical(c.planOnly) && isscalar(c.planOnly));
assert(islogical(c.makePlots) && isscalar(c.makePlots));
c.figureVisible = validatestring(c.figureVisible,{'on','off'});
c.derivativeDirections = string(c.derivativeDirections(:));
count = numel(c.familySets);
if isempty(c.derivativeDirections), c.derivativeDirections = repmat("none",count,1); end
assert(numel(c.derivativeDirections)==count);
base = bubble_options(input);
assert(base.n_level>0 && base.phi>0 && base.phi<pi, ...
    'The family sweep requires enhancing sources and a nondegenerate cap.');
c.maxMemoryGB = min(c.maxMemoryGB,base.maxSvdMemoryGB);
options = cell(count,1);
geometries = cell(count,1);
caseName = strings(count,1);
caseTag = strings(count,1);
columns = zeros(count,1);
naturalPoints = zeros(count,1);
eligible = false(count,1);
for k = 1:count
    opt = base;
    opt.refinementFamilies = c.familySets{k};
    opt.derivative_direction = c.derivativeDirections(k);
    opt.maxSvdMemoryGB = c.maxMemoryGB;
    options{k} = bubble_options(opt);
    g = bubble_geometry(options{k});
    geometries{k} = g;
    caseName(k) = strjoin(options{k}.refinementFamilies,'+');
    if opt.derivative_direction~="none"
        caseName(k) = caseName(k)+"+S_y("+opt.derivative_direction+")";
    end
    caseTag(k) = regexprep(caseName(k),'[^A-Za-z0-9]+','_');
    columns(k) = g.coefficientCount;
    naturalPoints(k) = g.collocationCount;
    estimate = bubble_svd_memory(3*naturalPoints(k),columns(k));
    eligible(k) = columns(k)<=c.maxSvdColumns && ...
        3*naturalPoints(k)<=c.maxSvdRows && ...
        naturalPoints(k)<=c.maxCollocationPoints && ...
        estimate.estimatedGB<=c.maxMemoryGB;
    assert(isequal(g.sourcePoints,geometries{1}.sourcePoints) && ...
        isequal(g.ringOffsets,geometries{1}.ringOffsets), ...
        'A family choice changed the source sites or ring offsets.');
end
assert(numel(unique(caseTag))==count,'Duplicate family cases.');
% Pick the densest affordable natural grid, then apply it to EVERY solved case.
% Reject larger cases explicitly rather than coarsening a grid to fit them.
sharedIndex = [];
if any(eligible)
    candidates = find(eligible);
    [~,index] = max(naturalPoints(candidates));
    sharedIndex = candidates(index);
    shared = geometries{sharedIndex};
    targetFields = {'points','normals','theta','cap','mobile','collocationCount', ...
        'capCount','mobileCount','baseCollocationCount','collocationRingIndex', ...
        'collocationRingCount','collocationRegion','ringOffsets','ringAzimuthCount', ...
        'edgeDistance','minimumEdgeDistance','minimumAddedEdgeDistance', ...
        'sourceRefinementCount'};
    for k = 1:count
        for j = 1:numel(targetFields)
            name = targetFields{j};
            geometries{k}.(name) = shared.(name);
        end
    end
end
fitPoints = zeros(count,1);
estimatedGB = zeros(count,1);
status = repmat("budget_limit",count,1);
for k = 1:count
    fitPoints(k) = geometries{k}.collocationCount;
    estimate = bubble_svd_memory(3*fitPoints(k),columns(k));
    estimatedGB(k) = estimate.estimatedGB;
    if eligible(k) && estimate.estimatedGB<=c.maxMemoryGB
        status(k) = "planned";
    end
end
summary = table(caseName,status,columns,naturalPoints,fitPoints,estimatedGB);
summary.retainedRank = nan(count,1);
summary.solveSeconds = nan(count,1);
summary.forceResidual = nan(count,1);
summary.relativeSpeedError = nan(count,1);
summary.fitMaxima = nan(count,4);
summary.coarseMaxima = nan(count,4);
summary.fineMaxima = nan(count,4);
summary.validationChange = nan(count,4);
summary.validationUnresolved = false(count,1);
summary.validatedMaxima = nan(count,4);
summary.overallMaximum = nan(count,1);
summary.ratioToStokeslets = nan(count,1);
result = struct('summary',summary,'controls',c,'caseOptions',{options}, ...
    'geometries',{geometries},'validation',[],'fineValidation',[], ...
    'profiles',{cell(count,1)},'lowestSampledCase',"",'bestResolvedCase',"");
disp(summary(:,1:6));
if c.planOnly, return; end
if strlength(string(c.outputDirectory))>0 && ~isfolder(c.outputDirectory)
    mkdir(c.outputDirectory);
end
if isempty(sharedIndex)
    save_summary(result);
    return
end
opt = options{sharedIndex};
offset = shared.ringOffsets(end);
validation = bubble_validation_grid(opt,shared.collocationCount,offset);
fineOpt = opt;
fineOpt.validationBaseRequest = max(opt.validationBaseRequest,2*shared.collocationCount);
fineOpt = bubble_options(fineOpt);
fineOffset = min(offset/2,validation.minimumEdgeDistance/opt.radius);
fineValidation = bubble_validation_grid(fineOpt,2*shared.collocationCount,fineOffset);
result.validation = validation;
result.fineValidation = fineValidation;
Fg = (4*pi*opt.radius^3/3)*opt.densityDifference*opt.gravity(:);
theoreticalUz = Fg(3)/stagnant_cap_exact_drag(opt.phi,0,opt.muOut,opt.radius);
for k = 1:count
    if summary.status(k)~="planned", continue; end
    opt = options{k};
    geometry = geometries{k};
    assert(isequal(geometry.points,shared.points));
    fprintf('\nFamily sweep %d/%d: %s\n',k,count,caseName(k));
    started = tic;
    solution = bubble_solve(opt,Fg,geometry);
    summary.solveSeconds(k) = toc(started);
    check = bubble_check_validation(solution,opt,validation,fineValidation);
    fitGrid = geometry;
    fitGrid.count = geometry.collocationCount;
    fitReport = bubble_validate(solution,fitGrid,opt);
    summary.status(k) = "completed";
    summary.retainedRank(k) = solution.retainedRank;
    summary.forceResidual(k) = solution.forceResidual;
    summary.relativeSpeedError(k) = abs(solution.U(3)/theoreticalUz-1);
    summary.fitMaxima(k,:) = fitReport.maxima;
    summary.coarseMaxima(k,:) = check.report.maxima;
    summary.fineMaxima(k,:) = check.fineReport.maxima;
    summary.validationChange(k,:) = check.relativeChange;
    summary.validationUnresolved(k) = check.unresolved;
    summary.validatedMaxima(k,:) = max(check.report.maxima,check.fineReport.maxima);
    summary.overallMaximum(k) = max(summary.validatedMaxima(k,:));
    result.profiles{k} = distance_profile(check.fineReport,fineValidation,opt);
    if strlength(string(c.outputDirectory))>0
        save(fullfile(c.outputDirectory,"oct29_bubble_families_"+caseTag(k)+"_data.mat"), ...
            'solution','opt','check','fitReport');
    end
    % Drop case-specific kernel/SVD data before solving the next family set.
    clear solution check fitReport
end
completed = find(summary.status=="completed");
baseline = find(caseName=="S" & summary.status=="completed",1);
if ~isempty(baseline)
    summary.ratioToStokeslets = summary.overallMaximum/summary.overallMaximum(baseline);
end
if ~isempty(completed)
    [~,index] = min(summary.overallMaximum(completed));
    result.lowestSampledCase = caseName(completed(index));
    resolved = completed(~summary.validationUnresolved(completed));
    if ~isempty(resolved)
        [~,index] = min(summary.overallMaximum(resolved));
        result.bestResolvedCase = caseName(resolved(index));
    end
    fprintf('\nLowest sampled worst-condition residual: %s.\n',result.lowestSampledCase);
    fprintf('Validation-stable winner (10%% grid-change test): %s\n',result.bestResolvedCase);
end
result.summary = summary;
disp(summary(:,{'caseName','status','columns','overallMaximum', ...
    'ratioToStokeslets','validationUnresolved'}));
save_summary(result);
if c.makePlots && ~isempty(completed)
    plot_comparison(result);
end
end

function save_summary(result)
% Preserve completed and budget-rejected outcomes without storing SVD factors.
directory = result.controls.outputDirectory;
if strlength(string(directory))==0, return; end
writetable(result.summary,fullfile(directory,'oct29_bubble_families_summary.csv'));
save(fullfile(directory,'oct29_bubble_families_summary.mat'),'result');
end

function profile = distance_profile(report,grid,opt)
% Bin worst residuals over all azimuths for comparable edge-distance curves.
distance = grid.edgeDistance/opt.radius;
edges = logspace(log10(min(distance))-1e-10,log10(max(distance))+1e-10,61);
bins = discretize(distance,edges);
values = [report.point.normal,report.point.normal, ...
    report.point.capTangent,report.point.mobileTraction]./report.normalization;
masks = [grid.cap,grid.mobile,grid.cap,grid.mobile];
profile.distanceOverA = sqrt(edges(1:end-1).*edges(2:end)).';
profile.maxima = nan(60,4);
for condition = 1:4
    for j = 1:60
        selected = masks(:,condition) & bins==j;
        if any(selected), profile.maxima(j,condition) = max(values(selected,condition)); end
    end
end
end

function plot_comparison(result)
% Overlay all completed bases using identical validation distances and scales.
c = result.controls;
f = figure('Color','w','Visible',c.figureVisible,'Position',[60 60 1350 850]);
layout = tiledlayout(f,2,2,'TileSpacing','compact');
titles = {'Cap normal velocity','Mobile normal velocity', ...
    'Cap tangential velocity','Mobile tangential traction'};
completed = find(result.summary.status=="completed");
for condition = 1:4
    ax = nexttile(layout);
    hold(ax,'on'); grid(ax,'on');
    set(ax,'XScale','log','YScale','log');
    for k = completed.'
        p = result.profiles{k};
        values = p.maxima(:,condition);
        values(values<1e-16) = 1e-16; % Preserve NaN gaps in empty distance bins.
        plot(ax,p.distanceOverA,values, ...
            'LineWidth',1.2,'DisplayName',result.summary.caseName(k));
    end
    yline(ax,1e-2,'k:','HandleVisibility','off');
    xlabel(ax,'Surface distance to edge s/a');
    ylabel(ax,'Normalized maximum per distance bin');
    title(ax,titles{condition});
    legend(ax,'Location','best','Interpreter','none');
end
sgtitle(layout,'Enhancing source families: shared fit and independent validation grids');
if strlength(string(c.outputDirectory))>0
    exportgraphics(f,fullfile(c.outputDirectory,'oct29_bubble_families_distance.png'), ...
        'Resolution',180);
end
end
