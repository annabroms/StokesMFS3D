function band = ring_band_geometry(centre,a,phi,depth,width,spacing,collocWidth,collocSpacing)
%RING_BAND_GEOMETRY Generate interior source and surface collocation bands.
% band = ring_band_geometry(centre,a,phi,depth,width,spacing,collocWidth,collocSpacing)
% centre is 1-by-3; the other inputs are scalars. Widths and target spacings
% are physical-surface arc lengths. Sources lie at radius a-depth; rings
% span phi +/- width/(2*a). Surface rings span collocWidth and avoid phi.
% band.sources is N-by-3; band.sourceRing is N-by-1, contiguous ring-major.
% band.collocation and band.normals are M-by-3, band.cap is M-by-1. Angles
% are colatitudes from +z, azimuths start at +x. Spacings are lengths, with
% projected source spacings measured at radius a. No strengths are normalized.
if nargin==0
    self_test();
    return
end
if ~isequal(size(centre),[1,3]) || any(~isfinite(centre)) || ...
        ~isscalar(a) || ~isfinite(a) || a<=0 || ...
        ~isscalar(phi) || ~isfinite(phi) || phi<=0 || phi>=pi || ...
        ~isscalar(depth) || ~isfinite(depth) || depth<=0 || depth>=a || ...
        ~isscalar(width) || ~isfinite(width) || width<=0 || ...
        ~isscalar(spacing) || ~isfinite(spacing) || spacing<=0 || ...
        ~isscalar(collocWidth) || ~isfinite(collocWidth) || collocWidth<=0 || ...
        ~isscalar(collocSpacing) || ~isfinite(collocSpacing) || collocSpacing<=0 || ...
        max(width,collocWidth)/(2*a)>=min(phi,pi-phi)
    error('ring_band_geometry:Geometry','Invalid band width, spacing, depth, sphere or cap angle.');
end
nr = max(1,ceil(width/spacing));
theta = phi-width/(2*a)+((1:nr)'-0.5)*width/(a*nr);
counts = max(3,round(2*pi*a*sin(theta)/spacing));
[unit,id] = sample_rings(theta,counts);
band.sources = centre+(a-depth)*unit;
if any(vecnorm(band.sources-centre,2,2)>=a)
    error('ring_band_geometry:SourceDepth', ...
        'The depth is too small to represent all sources strictly inside the sphere.');
end
band.sourceRing = id;
band.sourceTheta = theta;
band.sourceCounts = counts;
band.sourceMeridionalSpacing = width/nr;
band.sourceAzimuthSpacing = 2*pi*a*sin(theta)./counts;
band.depth = depth;

nc = max(1,ceil(collocWidth/(2*collocSpacing)));
half = collocWidth/(2*a);
offset = ((1:nc)'-0.5)*half/nc;
ctheta = [phi-flipud(offset);phi+offset];
ccounts = max(3,round(2*pi*a*sin(ctheta)/collocSpacing));
[normals,cid] = sample_rings(ctheta,ccounts);
band.collocation = centre+a*normals;
band.normals = normals;
band.collocRing = cid;
band.cap = normals(:,3)>cos(phi);
band.collocTheta = ctheta;
band.collocCounts = ccounts;
band.collocMeridionalSpacing = collocWidth/(2*nc);
band.collocAzimuthSpacing = 2*pi*a*sin(ctheta)./ccounts;
end

function [unit,id] = sample_rings(theta,counts)
%SAMPLE_RINGS Generate ordered azimuthal points for unequal ring counts.
N = sum(counts);
unit = zeros(N,3);
id = repelem((1:numel(theta))',counts);
start = 1;
for r = 1:numel(theta)
    rows = start:start+counts(r)-1;
    psi = (0:counts(r)-1)'*2*pi/counts(r);
    unit(rows,:) = [sin(theta(r))*cos(psi),sin(theta(r))*sin(psi), ...
        cos(theta(r))*ones(counts(r),1)];
    start = start+counts(r);
end
end

function self_test()
%SELF_TEST Verify geometry and show sources and surface-band collocation.
c = [0.2,-0.1,0.3]; a = 1; phi = pi/4;
b = ring_band_geometry(c,a,phi,0.2,0.6,0.13,0.8,0.18);
b2 = ring_band_geometry(c,a,phi,0.1,0.6,0.13,0.8,0.18);
assert(max(abs(vecnorm(b.sources-c,2,2)-0.8))<1e-12);
assert(max(abs(vecnorm(b.collocation-c,2,2)-a))<1e-12);
assert(max(abs(vecnorm(b.normals,2,2)-1))<1e-12);
assert(isequal(b.sourceCounts,b2.sourceCounts) && ...
    norm((b.sources-c)/0.8-(b2.sources-c)/0.9,'fro')<1e-12);
assert(sum(b.sourceCounts)==size(b.sources,1) && ...
    isequal(accumarray(b.sourceRing,1),b.sourceCounts));
assert(numel(unique(b.sourceCounts))>1 && any(b.cap) && any(~b.cap));
assert(all(abs(b.collocTheta-phi)>1e-12));
assert(all(b.sourceTheta>phi-0.3 & b.sourceTheta<phi+0.3));
try
    ring_band_geometry(c,a,phi,0,0.6,0.13,0.8,0.18);
    error('ring_band_geometry:SelfTest','Invalid depth accepted.');
catch err
    assert(strcmp(err.identifier,'ring_band_geometry:Geometry'));
end
try
    ring_band_geometry(c,a,phi,0.2,0.6,0.13,2,0.18);
    error('ring_band_geometry:SelfTest','Invalid surface-band width accepted.');
catch err
    assert(strcmp(err.identifier,'ring_band_geometry:Geometry'));
end
fig = figure('Color','w','Name','Source and collocation bands');
tl = tiledlayout(fig,1,2);
ax = nexttile(tl);
[sx,sy,sz] = sphere(32);
surf(ax,c(1)+a*sx,c(2)+a*sy,c(3)+a*sz, ...
    'FaceAlpha',0.08,'EdgeAlpha',0.08); hold(ax,'on');
scatter3(ax,b.sources(:,1),b.sources(:,2),b.sources(:,3),12,'b','filled');
scatter3(ax,b.collocation(:,1),b.collocation(:,2),b.collocation(:,3),12,'r','filled');
psi = linspace(0,2*pi,200);
plot3(ax,c(1)+a*sin(phi)*cos(psi),c(2)+a*sin(phi)*sin(psi), ...
    c(3)+a*cos(phi)*ones(size(psi)),'k','LineWidth',1.5);
axis(ax,'equal'); grid(ax,'on'); view(ax,35,22);
xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
title(ax,'Blue: interior sources; red: surface collocation');
ax = nexttile(tl); hold(ax,'on');
t = linspace(0,pi,300);
plot(ax,a*sin(t),a*cos(t),'k-',(a-0.2)*sin(t),(a-0.2)*cos(t),'b--');
plot(ax,(a-0.2)*sin(b.sourceTheta),(a-0.2)*cos(b.sourceTheta),'bo');
plot(ax,a*sin(b.collocTheta),a*cos(b.collocTheta),'r.');
plot(ax,a*sin(phi),a*cos(phi),'kx','MarkerSize',10,'LineWidth',2);
axis(ax,'equal'); grid(ax,'on'); xlabel(ax,'rho'); ylabel(ax,'z relative to centre');
title(ax,'Meridional section'); legend(ax,{'surface','source sphere','source rings', ...
    'collocation rings','edge'},'Location','best');
fprintf('ring_band_geometry self-test passed.\n');
end
