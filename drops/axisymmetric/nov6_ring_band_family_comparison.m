%% Compare mixed source families on an interior spherical band
% Run after setpath. All cases share source sites, surface collocation and
% independent check grids. Added surface band nodes increase their weight in
% the unweighted least-squares fit. See ring_band_geometry for conventions.

%% Physical, geometry and solver parameters
opt.a = 1; opt.mu = 1; opt.centre = [0 0 0]; opt.Fg = [0;0;-1];
opt.phi = pi/4; opt.ringMode = "swirl_free"; % "full" or "swirl_free"
opt.proxyRadius = 0.6*opt.a; opt.Nu = 32; opt.Nv = 20;
% The narrow band reduced SS_y residuals. Depth 0.07*a further improved
% the boundary maxima versus 0.035*a on a closer independent edge audit.
% Original setting: depth=0.05*a, sourceWidth=0.5*a, sourceSpacing=0.12*a.
depth = 0.07*opt.a; % use 0.035*a to recover the first narrow-band setting
opt.sourceWidth = 0.2*opt.a; opt.sourceSpacing = 0.035*opt.a;
opt.collocWidth = 0.7*opt.a; opt.collocSpacing = 0.1*opt.a;
opt.Mreq = 600; opt.collocationLevels = 8; opt.collocationSigma = 3;
opt.collocationAzimuths = 128; % [] uses the geometry-selected count
% Smaller angular offsets probe closer to the edge; 0.01 restores the old grid.
opt.validationMultiplier = 5; opt.validationOffset = 0.0004;
opt.auditFactor = 2; opt.svdTol = 1e-12; opt.beta = 1/(6*pi*opt.mu*opt.a);
% Empty is GL only. S, D, R, T_t, T_n, T_n2, Q_s and S_y are supported;
% S_y is the depth-scaled meridional source derivative of a Stokeslet.
familySets = {strings(1,0),"S","SD","ST_n","SS_y","ST_nS_y"};
mapCases = unique([1,min(2,numel(familySets)),numel(familySets)]); % edit as desired
if numel(familySets)<2 || ~isempty(familySets{1}) || ...
        ~isequal(ring_family_list(familySets{2}),"S")
    error('ring_band_family_comparison:Baselines','Keep GL-only and band-S as the first two cases.');
end

%% Shared geometry and independent check grids
data = ring_band_setup(opt,depth,familySets);
fprintf('Band depth %.6g, h/depth %.3g; interior spacing scale %.3g.\n', ...
    depth,opt.sourceSpacing/depth,(opt.a-depth)/opt.a);

%% Solve and validate on the same grids
results = cell(numel(familySets),1);
labels = strings(numel(familySets),1);
for j = 1:numel(familySets)
    results{j} = ring_family_case(data.Y0,data.id0,data.band.sources, ...
        data.band.sourceRing,opt.centre,opt.a,opt.mu,opt.Fg,familySets{j}, ...
        data.cg,data.op,data.val,data.audit,data.Uexact,opt.beta,opt.svdTol, ...
        opt.ringMode,opt.phi);
    labels(j) = results{j}.name;
end

%% Diagnostics and physical surface maps
ring_band_report(results,1:numel(results),labels,"Family combination", ...
    mapCases,norm(data.Uexact),opt.mu,opt.a);
for j = 3:numel(results)
    metrics = [results{j}.Uerror,results{j}.auditResidual];
    for k = 1:2
        baseline = [results{k}.Uerror,results{k}.auditResidual];
        fprintf('%s improves all four audit metrics over %s: %s\n', ...
            labels(j),labels(k),string(all(metrics<baseline)));
    end
end
