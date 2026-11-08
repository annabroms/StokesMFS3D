%% Compare mixed fundamental solutions on clustered source rings
% Run after setpath. All cases share the GL proxy, tapered collocation,
% clustered source sites and independent check grids. The background GL
% Stokeslets are always present; familySets specifies the additional blocks.
% Coefficients use rotating [radial; azimuthal; axial] components per ring.
% S, D and stresslet strength vectors are polar; their azimuthal components
% produce swirl for meridional source directions. A rotlet strength is axial:
% its azimuthal component produces meridional flow. T_n2 has azimuthal
% source direction, reversing these parity assignments for that stresslet.
% S has force units; R and T have force*length, D force*length^2.
% Q_s includes an a^2 kernel scale and has force-unit coefficients.
% S_y is a source-position Stokeslet derivative in the meridional direction,
% scaled by the source depth; it has force-unit coefficients.

%% Parameters
ringMode = "full";  % "full" or "swirl_free"
if ~isscalar(ringMode) || ~ismember(ringMode,["full","swirl_free"])
    error('ring_family_comparison:Mode','Use "full" or "swirl_free" for ringMode.');
end
phi = pi/4;
a = 1; mu = 1; centre = [0 0 0]; Fg = [0;0;-1];
proxyRadius = 0.6*a; Nu = 32; Nv = 20;
sourceLevels = 4; sourceStartDepth = 0.3*a; sourceSigma = 2;
% Azimuthal refinement: compare 128 with 256 while keeping the other settings
% fixed. Use 32 to recover the earlier source sampling. Ring unknowns are tied,
% so adding source points changes quadrature and cost, not the unknown count.
sourceAzimuths = 128;
Mreq = 600; collocationLevels = 8; collocationSigma = 3;
collocationAzimuths = 128; % [] restores the previous geometry-selected count
auditAzimuthFactor = 2;     % 1 restores the earlier source-count audit rings
validationMultiplier = 5; % stagnant_cap_validation_grid requires at least 5
validationOffset = 0.01;
svdTol = 1e-12; beta = 1/(6*pi*mu*a);
% Empty set is the GL-only baseline; S is the clustered-Stokeslet baseline.
% S Stokeslet, R rotlet, D potential dipole, T_t/T_n/T_n2 stresslets,
% Q_s axial quadrupole; S_y is the meridional source derivative.
% T is the legacy ray-oriented stresslet. All sources remain interior.
familySets = {strings(1,0),"S","SD","SR","ST_t","ST_tT_n","SQ_s","SS_y","ST_nS_y"};
momentWarningTol = 1e-6; % flag non-negligible force, torque or swirl residuals
if numel(familySets)<2 || ~isempty(familySets{1}) || ...
        ~isequal(ring_family_list(familySets{2}),"S")
    error('ring_family_comparison:Baselines', ...
        'The first two cases must be GL-only and clustered S.');
end

%% Shared geometry
[Y0,id0,nr] = ring_gl_sphere(proxyRadius,centre,Nu,Nv);
depth = stagnant_cap_cluster_offsets(sourceLevels,sourceStartDepth,"tapered",0.7,sourceSigma);
if any(depth>=a)
    error('ring_family_comparison:SourceDepth','Source depths must be less than the sphere radius.');
end
psi = (0:sourceAzimuths-1)'*2*pi/sourceAzimuths;
unit = [sin(phi)*cos(psi),sin(phi)*sin(psi),cos(phi)*ones(sourceAzimuths,1)];
Yc = zeros(sourceLevels*sourceAzimuths,3);
for j = 1:sourceLevels
    rows = (j-1)*sourceAzimuths+(1:sourceAzimuths);
    Yc(rows,:) = centre+(a-depth(j))*unit;
end
idc = repelem((1:sourceLevels)',sourceAzimuths);
for j = 1:size(Yc,1)
    if any(vecnorm(Y0-Yc(j,:),2,2)<1e-11*a)
        error('ring_family_comparison:CoincidentSources','GL and clustered sources coincide.');
    end
end
options = struct('radius',a,'centre',centre,'muOut',mu,'phi',phi,'sourceRequest',50, ...
    'collocationRatio',Mreq/size(stagnant_cap_design(50),1), ...
    'ringLevels',collocationLevels,'collocClustering',"tapered", ...
    'collocClusterSigma',collocationSigma,'ringAzimuthCount',collocationAzimuths);
cg = stagnant_cap_geometry(options);
M = cg.collocationCount;
maxUnknowns = 0;
for c = 2:numel(familySets)
    families = ring_family_list(familySets{c});
    if ringMode=="swirl_free"
        count = 2*numel(nr)+sourceLevels* ...
            (2*numel(families)-nnz(ismember(families,["R","T_n2"])));
    else
        count = 3*numel(nr)+3*sourceLevels*numel(families);
    end
    maxUnknowns = max(maxUnknowns,count);
end
if 3*M<maxUnknowns
    error('ring_family_comparison:Underdetermined', ...
        'The largest system has %d rows and %d coefficients; increase collocation.',3*M,maxUnknowns);
end
geom = struct('sourcePoints',Y0,'baseSourceCount',size(Y0,1), ...
    'coefficientCount',3*size(Y0,1),'collocationCount',M, ...
    'points',cg.points,'normals',cg.normals,'cap',cg.cap,'mobile',cg.mobile);
op = bubble_boundary_operators(geom,struct('n_level',0,'n_az_source',0,'muOut',mu));
gridOptions = struct('radius',a,'centre',centre,'muOut',mu,'phi',phi, ...
    'sourceRequest',50,'validationMultiplier',validationMultiplier);
val = stagnant_cap_validation_grid(gridOptions,M,validationOffset);
gridOptions.validationMultiplier = 2*validationMultiplier;
audit = stagnant_cap_validation_grid(gridOptions,M,validationOffset/2);
% The audit grid includes source-aligned and interlaced checks at each source
% depth and at its closest edge offset; these are not used in the solve.
auditAzimuths = auditAzimuthFactor*sourceAzimuths;
auditPsi = (0:auditAzimuths-1)'*2*pi/auditAzimuths;
auditOffsets = [depth(:).'/a,validationOffset/4];
for delta = auditOffsets
    for side = [-1,1]
        theta = phi+side*delta;
        if theta<=0 || theta>=pi, continue; end
        for phase = [0,0.5]
            angle = auditPsi+phase*2*pi/auditAzimuths;
            n = [sin(theta)*cos(angle),sin(theta)*sin(angle),cos(theta)*ones(size(angle))];
            audit.points = [audit.points;centre+a*n];
            audit.normals = [audit.normals;n];
            audit.cap = [audit.cap;repmat(theta<phi,auditAzimuths,1)];
        end
    end
end
audit.mobile = ~audit.cap;
audit.theta = acos(max(-1,min(1,audit.normals(:,3))));
audit.edgeDistance = a*abs(audit.theta-phi);
audit.nearEdge = audit.edgeDistance<=a*cg.bandAngle;
audit.minimumEdgeDistance = min(audit.edgeDistance);
audit.count = size(audit.points,1);
fprintf('Mode %s; GL: %d rings, %d sources; clustered: %d rings, %d sources.\n', ...
    ringMode,numel(nr),size(Y0,1),sourceLevels,size(Yc,1));
fprintf('Source depths:'); fprintf(' %.3g',depth); fprintf('\n');
fprintf('Collocation %d, check %d, audit %d; closest source depth %.3g.\n', ...
    M,val.count,audit.count,min(depth));
fprintf('Azimuths per ring: source %d, collocation %d, extra audit %d per phase.\n', ...
    sourceAzimuths,cg.ringAzimuthCount,auditAzimuths);
arcSpacing = 2*pi*a*sin(phi)/sourceAzimuths;
if arcSpacing>min(depth)
    warning('ring_family_comparison:AzimuthResolution', ...
        'Source azimuth spacing %.3g exceeds minimum source depth %.3g; audit for narrow spikes.', ...
        arcSpacing,min(depth));
end

%% Solve and validate every family set
zeta = stagnant_cap_exact_drag(phi,0,mu,a);
Uexact = Fg/zeta; Uref = norm(Uexact); tractionScale = mu*Uref/a;
result = cell(numel(familySets),1);
for c = 1:numel(familySets)
    result{c} = ring_family_case(Y0,id0,Yc,idc,centre,a,mu,Fg, ...
        familySets{c},cg,op,val,audit,Uexact,beta,svdTol,ringMode);
end

%% Report and plots
names = strings(numel(result),1);
metrics = zeros(numel(result),4); auditMetrics = metrics;
fprintf('\nCase                       U error    cap vel    mobile n   traction    force      torque     swirl\n');
for c = 1:numel(result)
    r = result{c}; names(c) = r.name;
    metrics(c,:) = [r.Uerror,r.residual];
    auditMetrics(c,:) = [r.Uerror,r.auditResidual];
    fprintf('%-26s %.2e   %.2e   %.2e   %.2e   %.2e   %.2e   %.2e\n', ...
        r.name,metrics(c,:),r.forceResidual,r.torqueResidual,r.swirl);
    fprintf('  audit boundary maxima: %.2e %.2e %.2e; rank %d/%d\n', ...
        r.auditResidual,r.rank,r.unknowns);
    if any(r.auditResidual>1.5*max(r.residual,eps))
        fprintf('  WARNING: denser audit increases at least one maximum by more than 50%%.\n');
    end
    fprintf('  regional maxima [near edge, far edge, polar], rows cap/normal/traction:\n');
    disp(r.regional);
    if max([r.forceResidual,r.torqueResidual,r.swirl])>momentWarningTol
        fprintf('  WARNING: force, torque or swirl exceeds %.1e; assess physical reliability.\n',momentWarningTol);
    end
    for b = 1:numel(r.familyNames)
        fprintf('  %s max physical coefficient / scale = %.3e (scale %.3e)\n', ...
            r.familyNames(b),r.coefficientMax(b),r.coefficientScale(b));
    end
end
for c = 3:numel(result)
    for base = 1:2
        improves = all(auditMetrics(c,:)<auditMetrics(base,:));
        fprintf('%s improves all four metrics over %s: %s\n', ...
            names(c),names(base),string(improves));
    end
end
figure('Color','w','Position',[40 80 1400 500]);
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl);
semilogy(ax,1:numel(result),max(auditMetrics,eps),'o-');
grid(ax,'on'); xticks(ax,1:numel(result)); xticklabels(ax,names); xtickangle(ax,45);
legend(ax,{'U error','cap velocity','mobile normal','mobile traction'},'Location','best');
ylabel(ax,'relative error / normalized maximum'); title(ax,'Denser independent audit');
ax = nexttile(tl);
semilogy(ax,1:numel(result),max(cellfun(@(r) r.forceResidual,result),eps),'o-', ...
    1:numel(result),max(cellfun(@(r) r.torqueResidual,result),eps),'s-', ...
    1:numel(result),max(cellfun(@(r) r.swirl,result),eps),'^-');
grid(ax,'on'); xticks(ax,1:numel(result)); xticklabels(ax,names); xtickangle(ax,45);
legend(ax,{'force / |F_g|','torque / (a|F_g|)','swirl / U_{ex}'},'Location','best');
title(ax,'Moment and swirl diagnostics');
ax = nexttile(tl); hold(ax,'on');
allFamilies = strings(1,0);
for c = 1:numel(result)
    allFamilies = union(allFamilies,result{c}.familyNames,'stable');
end
for j = 1:numel(allFamilies)
    values = nan(numel(result),1);
    for c = 1:numel(result)
        index = find(result{c}.familyNames==allFamilies(j),1);
        if ~isempty(index), values(c) = result{c}.coefficientMax(index); end
    end
    values(values==0) = eps;
    semilogy(ax,1:numel(result),values,'o','DisplayName',allFamilies(j));
end
grid(ax,'on'); xticks(ax,1:numel(result));
xticklabels(ax,names); xtickangle(ax,45);
ylabel(ax,'max per-node coefficient / family scale');
title(ax,'Family-wise coefficient magnitude');
legend(ax,'Location','best');
chosen = find(all(auditMetrics(3:end,:)<auditMetrics(2,:),2),1)+2;
if isempty(chosen), chosen = min(3,numel(result)); end
for c = unique([1,2,chosen])
    r = result{c};
    ring_plot_residuals(r.mapPoints,r.mapResiduals, ...
        {char(r.name+" cap |u-U|/Uref"),char(r.name+" mobile |(u-U).n|/Uref"), ...
        char(r.name+" mobile tangential traction/(mu Uref/a)")});
end
fprintf('Uref=%.6g; traction scale=%.6g; plotting uses independent audit nodes.\n',Uref,tractionScale);
