function data = ring_band_setup(opt,depths,familySets)
%RING_BAND_SETUP Prepare common band, collocation and independent audit grids.
% opt contains the parameters used by the band experiments; depths is a
% positive vector and familySets is a cell array of ring_family_list tokens.
% data contains GL and band coordinates, ring IDs, shared geometry/operators,
% check and audit grids, and the exact axial translation. Widths and spacings
% are physical-surface lengths; normals point outward and z is the polar axis.
% Additional surface rows are unweighted and increase least-squares emphasis
% on the band. Node vectors are interleaved [x1;y1;z1;...].
if nargin==0
    self_test();
    return
end
a = opt.a; centre = opt.centre; phi = opt.phi;
if isempty(depths) || any(~isfinite(depths)) || any(depths<=0 | depths>=a)
    error('ring_band_setup:Depth','Every band depth must lie strictly between zero and a.');
end
if ~isscalar(opt.ringMode) || ~ismember(opt.ringMode,["full","swirl_free"]) || ...
        ~isscalar(opt.auditFactor) || ~isfinite(opt.auditFactor) || ...
        opt.auditFactor<1 || opt.auditFactor~=round(opt.auditFactor)
    error('ring_band_setup:Options','Use a supported ring mode and a positive integer audit factor.');
end
[Y0,id0,nr] = ring_gl_sphere(opt.proxyRadius,centre,opt.Nu,opt.Nv);
band = ring_band_geometry(centre,a,phi,min(depths),opt.sourceWidth, ...
    opt.sourceSpacing,opt.collocWidth,opt.collocSpacing);
for d = depths(:).'
    Y = centre+(a-d)/(a-min(depths))*(band.sources-centre);
    if any(vecnorm(Y-centre,2,2)>=a)
        error('ring_band_setup:Depth','A depth is too small to keep sources numerically interior.');
    end
    for j = 1:size(Y,1)
        if any(vecnorm(Y0-Y(j,:),2,2)<1e-11*a)
            error('ring_band_setup:CoincidentSources','Band and GL sources coincide.');
        end
    end
end
options = struct('radius',a,'centre',centre,'muOut',opt.mu,'phi',phi,'sourceRequest',50, ...
    'collocationRatio',opt.Mreq/size(stagnant_cap_design(50),1), ...
    'ringLevels',opt.collocationLevels,'collocClustering',"tapered", ...
    'collocClusterSigma',opt.collocationSigma,'ringAzimuthCount',opt.collocationAzimuths);
cg = stagnant_cap_geometry(options);
% Preserve old rows and discard only proposed band nodes within 1e-11*a.
tol = 1e-11*a;
[points,accepted] = append_distinct(cg.points,band.collocation,tol);
oldCount = cg.collocationCount;
cg.points = points;
cg.normals = (points-centre)/a;
cg.theta = acos(max(-1,min(1,cg.normals(:,3))));
cg.cap = cg.theta<phi;
cg.mobile = ~cg.cap;
cg.edgeDistance = a*abs(cg.theta-phi);
cg.collocationCount = size(points,1);
cg.baseCollocation = [cg.baseCollocation;false(nnz(accepted),1)];
cg.ringCollocation = [cg.ringCollocation;false(nnz(accepted),1)];
cg.bandCollocation = [false(oldCount,1);true(nnz(accepted),1)];
if any(cg.edgeDistance<=tol)
    error('ring_band_setup:EdgeCollocation','Collocation must not lie on the cap edge.');
end
R = numel(band.sourceCounts);
maxUnknowns = 0;
for c = 1:numel(familySets)
    families = ring_family_list(familySets{c});
    if opt.ringMode=="swirl_free"
        count = 2*numel(nr)+R*(2*numel(families)-nnz(ismember(families,["R","T_n2"])));
    else
        count = 3*numel(nr)+3*R*numel(families);
    end
    maxUnknowns = max(maxUnknowns,count);
end
if 3*cg.collocationCount<maxUnknowns
    error('ring_band_setup:Underdetermined','Increase collocation: %d rows for %d coefficients.', ...
        3*cg.collocationCount,maxUnknowns);
end
geom = struct('sourcePoints',Y0,'baseSourceCount',size(Y0,1), ...
    'coefficientCount',3*size(Y0,1),'collocationCount',cg.collocationCount, ...
    'points',cg.points,'normals',cg.normals,'cap',cg.cap,'mobile',cg.mobile);
op = bubble_boundary_operators(geom,struct('n_level',0,'n_az_source',0,'muOut',opt.mu));
gridOpt = struct('radius',a,'centre',centre,'muOut',opt.mu,'phi',phi, ...
    'sourceRequest',50,'validationMultiplier',opt.validationMultiplier);
val = stagnant_cap_validation_grid(gridOpt,cg.collocationCount,opt.validationOffset);
gridOpt.validationMultiplier = 2*opt.validationMultiplier;
audit = stagnant_cap_validation_grid(gridOpt,cg.collocationCount,opt.validationOffset/2);
% Resolve each source projection and the surrounding footprints using
% aligned/interlaced azimuths and offsets at the smallest source depth.
delta = min(depths)/(2*a);
for r = 1:R
    for theta = band.sourceTheta(r)+[-delta,0,delta]
        if theta<=0 || theta>=pi || abs(theta-phi)<=1e-12, continue; end
        audit = append_ring(audit,centre,a,phi,theta,opt.auditFactor*band.sourceCounts(r));
    end
end
% Interlaced checks cover the independently specified collocation band.
for r = 1:numel(band.collocTheta)
    theta = band.collocTheta(r);
    audit = append_ring(audit,centre,a,phi,theta,opt.auditFactor*band.collocCounts(r));
end
[val,removedVal] = exclude_collocation(val,cg.points,tol,a,phi,cg.bandAngle);
[audit,removedAudit] = exclude_collocation(audit,cg.points,tol,a,phi,cg.bandAngle);
audit.mobile = ~audit.cap;
audit.theta = acos(max(-1,min(1,audit.normals(:,3))));
audit.edgeDistance = a*abs(audit.theta-phi);
audit.nearEdge = audit.edgeDistance<=a*cg.bandAngle;
audit.minimumEdgeDistance = min(audit.edgeDistance);
audit.count = size(audit.points,1);
zeta = stagnant_cap_exact_drag(phi,0,opt.mu,a);
data = struct('Y0',Y0,'id0',id0,'band',band,'cg',cg,'op',op, ...
    'val',val,'audit',audit,'Uexact',opt.Fg/zeta);
fprintf('GL %d rings / %d sources; band %d rings / %d sources.\n', ...
    numel(nr),size(Y0,1),R,size(band.sources,1));
fprintf('Source projected spacing: meridional %.3g, azimuth %.3g to %.3g.\n', ...
    band.sourceMeridionalSpacing,min(band.sourceAzimuthSpacing),max(band.sourceAzimuthSpacing));
fprintf('Collocation: base %d, edge %d, band %d (%d duplicates omitted), total %d.\n', ...
    nnz(cg.baseCollocation),nnz(cg.ringCollocation),nnz(cg.bandCollocation), ...
    nnz(~accepted),cg.collocationCount);
for mask = {cg.baseCollocation,cg.ringCollocation,cg.bandCollocation}
    fprintf('  cap/mobile: %d/%d\n',nnz(mask{1} & cg.cap),nnz(mask{1} & cg.mobile));
end
fprintf('Collocation band width %.3g; spacing meridional %.3g, azimuth %.3g to %.3g.\n', ...
    opt.collocWidth,band.collocMeridionalSpacing, ...
    min(band.collocAzimuthSpacing),max(band.collocAzimuthSpacing));
fprintf('Independent check %d, common audit %d; removed collocation overlaps %d/%d.\n', ...
    val.count,audit.count,removedVal,removedAudit);
fprintf('Minimum source depth %.3g.\n',min(depths));
auditSpacing = max([band.sourceMeridionalSpacing,band.collocMeridionalSpacing, ...
    max(band.sourceAzimuthSpacing)/opt.auditFactor, ...
    max(band.collocAzimuthSpacing)/opt.auditFactor]);
if max([band.sourceMeridionalSpacing;band.sourceAzimuthSpacing])>min(depths) || ...
        auditSpacing>min(depths)
    warning('ring_band_setup:Resolution', ...
        'Source or audit spacing exceeds the minimum depth; narrow residual spikes may be missed.');
end
end

function [points,accepted] = append_distinct(points,added,tol)
%APPEND_DISTINCT Preserve existing points while rejecting overlapping additions.
accepted = false(size(added,1),1);
for j = 1:size(added,1)
    if all(vecnorm(points-added(j,:),2,2)>tol)
        accepted(j) = true;
        points(end+1,:) = added(j,:); %#ok<AGROW>
    end
end
end

function [grid,removed] = exclude_collocation(grid,collocation,tol,a,phi,bandAngle)
%EXCLUDE_COLLOCATION Ensure validation nodes are independent of solve rows.
overlap = ismembertol(grid.points,collocation,tol,'ByRows',true,'DataScale',1);
removed = nnz(overlap);
grid.baseCount = grid.baseCount-nnz(overlap(1:grid.baseCount));
grid.points(overlap,:) = [];
grid.normals(overlap,:) = [];
grid.cap(overlap,:) = [];
grid.mobile = ~grid.cap;
grid.theta = acos(max(-1,min(1,grid.normals(:,3))));
grid.edgeDistance = a*abs(grid.theta-phi);
grid.nearEdge = grid.edgeDistance<=a*bandAngle;
grid.count = size(grid.points,1);
grid.minimumEdgeDistance = min(grid.edgeDistance);
end

function grid = append_ring(grid,centre,a,phi,theta,nAz)
%APPEND_RING Add source-aligned and interlaced physical-surface audit points.
psi = (0:nAz-1)'*2*pi/nAz;
for phase = [0,0.5]
    angle = psi+phase*2*pi/nAz;
    n = [sin(theta)*cos(angle),sin(theta)*sin(angle),cos(theta)*ones(nAz,1)];
    grid.points = [grid.points;centre+a*n];
    grid.normals = [grid.normals;n];
    grid.cap = [grid.cap;repmat(theta<phi,nAz,1)];
end
end

function self_test()
%SELF_TEST Check combined-grid bookkeeping and independent masks.
o = struct('a',1,'centre',[0 0 0],'phi',pi/4,'mu',1,'Fg',[0;0;-1], ...
    'proxyRadius',0.45,'Nu',6,'Nv',4,'sourceWidth',0.4,'sourceSpacing',0.25, ...
    'collocWidth',0.5,'collocSpacing',0.3,'Mreq',60,'collocationLevels',2, ...
    'collocationSigma',2,'collocationAzimuths',12,'ringMode',"swirl_free", ...
    'validationMultiplier',5,'validationOffset',0.01,'auditFactor',2);
d = ring_band_setup(o,[0.15,0.2],{strings(1,0),"SS_y"});
g = d.cg;
assert(sum(g.baseCollocation | g.ringCollocation | g.bandCollocation)==g.collocationCount);
assert(~any(g.baseCollocation & g.ringCollocation | g.bandCollocation & g.ringCollocation));
assert(nnz(g.bandCollocation)>0 && all(g.edgeDistance>1e-11));
assert(any(g.cap & g.bandCollocation) && any(g.mobile & g.bandCollocation));
assert(all(d.audit.edgeDistance>0) && d.audit.count>g.collocationCount);
assert(size(d.op.Nm,1)==nnz(g.mobile));
[combined,keep] = append_distinct([0 0 0],[0 0 0;1 0 0;1+1e-12 0 0],1e-11);
assert(isequal(keep,[false;true;false]) && size(combined,1)==2);
assert(~any(ismembertol(d.val.points,g.points,1e-11,'ByRows',true,'DataScale',1)));
assert(~any(ismembertol(d.audit.points,g.points,1e-11,'ByRows',true,'DataScale',1)));
fprintf('ring_band_setup self-test passed.\n');
end
