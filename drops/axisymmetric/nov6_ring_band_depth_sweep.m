%% Sweep interior source-band depth for one fixed family combination
% Run after setpath. Source angles and ring populations, collocation and both
% independent check grids are fixed throughout the sweep. Added surface band
% nodes increase their weight in the unweighted least-squares fit.

%% Physical, geometry and solver parameters
opt.a = 1; opt.mu = 1; opt.centre = [0 0 0]; opt.Fg = [0;0;-1];
opt.phi = pi/4; opt.ringMode = "swirl_free"; % "full" or "swirl_free"
opt.proxyRadius = 0.6*opt.a; opt.Nu = 32; opt.Nv = 20;
% Tested close bands; former depth list was [0.08 0.12 0.18 0.25].
depths = opt.a*[0.035 0.05 0.07 0.08 0.12 0.18 0.25];
% Previous sourceWidth=0.5*a and sourceSpacing=0.12*a.
opt.sourceWidth = 0.2*opt.a; opt.sourceSpacing = 0.035*opt.a;
opt.collocWidth = 0.7*opt.a; opt.collocSpacing = 0.1*opt.a;
opt.Mreq = 600; opt.collocationLevels = 8; opt.collocationSigma = 3;
opt.collocationAzimuths = 128;
% Probe closer to the edge; 0.01 restores the old validation offset.
opt.validationMultiplier = 5; opt.validationOffset = 0.0004;
opt.auditFactor = 2; opt.svdTol = 1e-12; opt.beta = 1/(6*pi*opt.mu*opt.a);
families = "SS_y"; % fixed throughout; S_y is a meridional source derivative
mapCases = unique([1,ceil(numel(depths)/2),numel(depths)]); % edit as desired

%% Shared angular geometry and independent check grids
data = ring_band_setup(opt,depths,{families});

%% Change only radial source distance, then solve on common grids
results = cell(numel(depths),1);
labels = strings(numel(depths),1);
for j = 1:numel(depths)
    depth = depths(j);
    Y = opt.centre+(opt.a-depth)/(opt.a-data.band.depth)* ...
        (data.band.sources-opt.centre);
    fprintf('\nDepth %.6g, h/depth %.3g; interior spacing scale %.3g.\n', ...
        depth,opt.sourceSpacing/depth,(opt.a-depth)/opt.a);
    results{j} = ring_family_case(data.Y0,data.id0,Y,data.band.sourceRing, ...
        opt.centre,opt.a,opt.mu,opt.Fg,families,data.cg,data.op,data.val, ...
        data.audit,data.Uexact,opt.beta,opt.svdTol,opt.ringMode,opt.phi);
    labels(j) = "depth/a = "+string(depth/opt.a);
end

%% Diagnostics and physical surface maps
ring_band_report(results,depths/opt.a,labels,"Depth / a",mapCases, ...
    norm(data.Uexact),opt.mu,opt.a);
