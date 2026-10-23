function results = compare_drop_image_source_families(config)
%COMPARE_DROP_IMAGE_SOURCE_FAMILIES Compare image bases at gap delta=0.1.
% Run with no input for the pressure-study geometry at requested proxy counts
% 700 and 1000. This is a substantial experiment; a small config can be used
% as a smoke check. A difference above 1e-3 measures basis sensitivity.
% imageFamilies tokens: S Stokeslet, R rotlet, D potential dipole, T legacy
% stresslet, T_t line-parallel, T_n first transverse, T_n2 second transverse
% stresslet, Q_s trace-free axial Stokeslet quadrupole. T excludes T_*.
% The four ordered selections below keep the original SDRT_nT_t basis intact.

if nargin==0, config = struct(); end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(root));
defaults = struct('requestedN',[700 1000],'validationCount',4000, ...
    'validationGapCounts',[4001 1001],'offSurfaceCount',50, ...
    'capPointCounts',[800 200],'N_image',30, ...
    'directions',[0 0 1;0 0 1], ...
    'families',{{"SDRT_nT_t","SDRT_nT_tT_n2", ...
    "SDRT_nT_tQ_s","SDRT_nT_tT_n2Q_s"}});
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(config,names{k}), config.(names{k}) = defaults.(names{k}); end
end
q = [-1.05 0 0;1.05 0 0]; a = ones(2,1); S1 = 2*ones(2,1);
e = config.directions;
assert(isequal(size(e),[2 3]) && all(vecnorm(e,2,2)>0));
e = e./vecnorm(e,2,2);
families = string(config.families);
assert(numel(families)==4 && families(1)=="SDRT_nT_t");
for k = 1:4, image_family_list(families(k)); end
fieldScale = max(abs(S1));
velocityScale = max(abs(S1).*a);
opt = struct('N',config.requestedN(1),'proxyRadius',0.63, ...
    'outerProxyRatio',2,'collocationRatio',2,'useImages',true, ...
    'N_image',config.N_image,'imageFamilies',families(1), ... % See token legend above.
    'mirror_ext',true,'nearThreshold',0.15, ...
    'outerImageOffsetFraction',0.001,'capMode',"exponential", ...
    'capAngles',[pi/5 pi/60],'capPointCounts',config.capPointCounts, ...
    'capExponent',[3 3],'capLayers',[18 18], ...
    'useRowScaling',true,'rowWeightMode',"local", ...
    'interiorViscosity',10,'completionScale',1/(8*pi), ...
    'velocityScale',velocityScale,'tractionScale',fieldScale, ...
    'tractionRowWeight',1,'scale_columns',false,'svdCutoff',1e-12, ...
    'gmresTolerance',1e-10,'maxIterations',350, ...
    'batchSize',128,'solverMode',"oneBody",'plot_grid',false);
opt = drop_image_options(opt,2,a);

results = struct('NRequested',{},'NActual',{},'family',{}, ...
    'region',{},'location',{},'velocityDifference',{}, ...
    'pressureDifference',{},'tractionDifference',{}, ...
    'continuityResidual',{},'normalResidual',{}, ...
    'tangentialResidual',{},'retainedRanks',{}, ...
    'systemSize',{},'svdCutoff',{},'gmresIterations',{}, ...
    'gmresResidual',{},'innerProxyMax',{},'outerProxyMax',{}, ...
    'innerImageMax',{},'outerImageMax',{},'pressureGauge',{});
row = 0;
for ni = 1:numel(config.requestedN)
    opt.N = config.requestedN(ni);
    referenceGrid = drop_image_discretization(q,a,opt);
    [surface,surfaceNormals] = drop_image_validation_grid(referenceGrid, ...
        config.validationCount,config.validationGapCounts,[pi/5 pi/60]);
    [off,offNormals] = independent_off_surface_points( ...
        q,a,config.offSurfaceCount);
    reference = [];
    for k = 1:numel(families)
        opt.imageFamilies = families(k); % See token legend above.
        fprintf('\nN requested %d, family %s, directions [%g %g %g]; [%g %g %g]\n', ...
            opt.N,families(k),e.');
        [U,lambda,gamma,grid,info] = solve_drops_with_images(q,a,S1,e,opt);
        for body = 1:2
            assert(isequal(grid.body(body).collocation, ...
                referenceGrid.body(body).collocation));
        end
        maxima = drop_image_validate_surface_conditions( ...
            grid,lambda,gamma,U,S1,e,info.options,surface,surfaceNormals);
        maxima = max(maxima,[],1)./[velocityScale velocityScale fieldScale];
        data = evaluate_fields(grid,lambda,gamma,info.options, ...
            surface,surfaceNormals,off,offNormals);
        [innerProxy,outerProxy,innerImages,outerImages] = ...
            drop_image_max_source_strengths(grid,lambda,gamma,families(k));
        if k==1, reference = data; end
        for region = 0:2
            gauge = 0;
            if region>0
                gauge = mean(data(region+1).off.p- ...
                    reference(region+1).off.p);
            end
            for location = ["surface","off"]
                candidate = data(region+1).(location);
                baseline = reference(region+1).(location);
                dp = candidate.p-baseline.p-gauge;
                dt = candidate.t-baseline.t+gauge*candidate.n;
                row = row+1;
                results(row).NRequested = opt.N;
                results(row).NActual = grid.N;
                results(row).family = families(k);
                results(row).region = region;
                results(row).location = location;
                results(row).velocityDifference = max( ...
                    vecnorm(candidate.u-baseline.u,2,2))/velocityScale;
                results(row).pressureDifference = max(abs(dp))/fieldScale;
                results(row).tractionDifference = max(vecnorm(dt,2,2))/fieldScale;
                results(row).continuityResidual = maxima(1);
                results(row).normalResidual = maxima(2);
                results(row).tangentialResidual = maxima(3);
                results(row).retainedRanks = info.retainedRanks;
                results(row).systemSize = info.systemSize;
                results(row).svdCutoff = opt.svdCutoff;
                results(row).gmresIterations = info.iterations;
                results(row).gmresResidual = info.originalResidual;
                results(row).innerProxyMax = innerProxy;
                results(row).outerProxyMax = outerProxy;
                results(row).innerImageMax = innerImages;
                results(row).outerImageMax = outerImages;
                results(row).pressureGauge = gauge;
            end
        end
        print_case(results(row-5:row));
    end
end
end

function [off,normals] = independent_off_surface_points(q,a,count)
% Place common gap targets in the exterior and radial targets in each drop.
[v,~] = get_sphdesign(count);
v = v./vecnorm(v,2,2);
off = cell(3,1); normals = cell(3,1);
off{1} = [zeros(size(v,1),1),0.55*v(:,2:3)];
normals{1} = v;
for body = 1:2
    off{body+1} = q(body,:)+0.7*a(body)*v;
    normals{body+1} = v;
end
end

function data = evaluate_fields(grid,lambda,gamma,opt,surface,surfaceNormals,off,offNormals)
% Evaluate all common points and three off-surface traction normals.
data = repmat(struct('surface',[],'off',[]),3,1);
for region = 0:2
    if region==0
        x = vertcat(surface{:}); n = vertcat(surfaceNormals{:});
    else
        x = surface{region}; n = surfaceNormals{region};
    end
    [u,p,t] = drop_image_field(x,n,region,grid,lambda,gamma,opt);
    data(region+1).surface = struct('u',u,'p',p,'t',t,'n',n);
    x = off{region+1}; n = offNormals{region+1};
    [u,p,t] = drop_image_field(x,n,region,grid,lambda,gamma,opt);
    % Use three independent normals so the off-surface result tests all
    % columns of the Cauchy stress, not only radial traction.
    for axis = 1:3
        nn = repmat(double((1:3)==axis),size(x,1),1);
        [~,~,extra] = drop_image_field(x,nn,region,grid,lambda,gamma,opt);
        t = [t;extra]; %#ok<AGROW>
        n = [n;nn]; %#ok<AGROW>
    end
    data(region+1).off = struct('u',u,'p',p, ...
        't',t,'n',n);
end
end

function print_case(rows)
% Print scaled maxima and solver diagnostics for one source selection.
r = rows(1);
fprintf('  actual N %d, ranks %s, SVD cutoff %.0e, GMRES %d, residual %.2e\n', ...
    r.NActual,mat2str(r.retainedRanks),r.svdCutoff, ...
    r.gmresIterations,r.gmresResidual);
fprintf('  BC [continuity normal tangent] = [%.2e %.2e %.2e]\n', ...
    r.continuityResidual,r.normalResidual,r.tangentialResidual);
fprintf('  max strengths inner/outer proxy %.2e %.2e; images %s / %s\n', ...
    r.innerProxyMax,r.outerProxyMax, ...
    mat2str(r.innerImageMax,3),mat2str(r.outerImageMax,3));
for j = 1:numel(rows)
    r = rows(j);
    fprintf('  region %d %s: du %.2e, dp %.2e, dt %.2e, gauge %.2e\n', ...
        r.region,r.location,r.velocityDifference, ...
        r.pressureDifference,r.tractionDifference,r.pressureGauge);
end
end
