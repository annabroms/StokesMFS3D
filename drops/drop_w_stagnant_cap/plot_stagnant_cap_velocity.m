function [surfaceFigure,profileFigure,diagnostic] = ...
    plot_stagnant_cap_velocity(solution,opt,visibility)
%PLOT_STAGNANT_CAP_VELOCITY Plot drop-relative meridional surface velocity.
if nargin == 0
    self_test();
    return
end
if nargin<3, visibility = 'on'; end
opt = stagnant_cap_options(opt);
if isfield(solution,'exteriorBasis')
    exteriorBasis = solution.exteriorBasis;
else
    exteriorBasis = struct('family',"S", ...
        'points',solution.grid.innerSources, ...
        'normalDirections',zeros(0,3),'scale',1);
end
assert(opt.phi>0 && opt.phi<pi, ...
    'Edge profile requires a nondegenerate stagnant cap.');
a = opt.radius; c = opt.centre;
theta = ((0:opt.plotThetaCount-1)+0.5)*pi/opt.plotThetaCount;
psi = 2*pi*(0:opt.plotAzimuthCount-1)/opt.plotAzimuthCount;
[thetaGrid,psiGrid] = ndgrid(theta,psi);
unit = [sin(thetaGrid(:)).*cos(psiGrid(:)), ...
    sin(thetaGrid(:)).*sin(psiGrid(:)),cos(thetaGrid(:))];
eTheta = [cos(thetaGrid(:)).*cos(psiGrid(:)), ...
    cos(thetaGrid(:)).*sin(psiGrid(:)),-sin(thetaGrid(:))];
points = c+a*unit;
relative = zeros(size(points));
for first = 1:opt.batchSize:size(points,1)
    rows = first:min(first+opt.batchSize-1,size(points,1));
    u = stokes_mfs_velocity_matrix( ...
        exteriorBasis,points(rows,:))*(solution.lambda/opt.muOut);
    relative(rows,:) = reshape(u,3,[]).'-solution.U(:).';
end
uTheta = sum(relative.*eTheta,2);
plot_surface_scalar(points,size(points,1),1,uTheta, ...
    'Meridional velocity relative to translating drop','parula');
surfaceFigure = gcf;
set(surfaceFigure,'Visible',visibility,'Name','Stagnant-cap surface velocity');
hold(gca,'on');
edgePsi = linspace(0,2*pi,400);
plot3(gca,c(1)+a*sin(opt.phi)*cos(edgePsi), ...
    c(2)+a*sin(opt.phi)*sin(edgePsi), ...
    c(3)+a*cos(opt.phi)*ones(size(edgePsi)), ...
    'k-','LineWidth',1.6);
cb = findall(surfaceFigure,'Type','ColorBar');
if ~isempty(cb), cb(1).Label.String = 'u_\theta-U\cdot e_\theta'; end

baseTheta = ((0:159)+0.5)*pi/160;
if opt.phi>0 && opt.phi<pi && ~isempty(solution.grid.ringOffsets)
    smallest = solution.grid.ringOffsets(end)/4;
    largest = min(opt.edgeBandAngle,pi-opt.phi);
    offsets = logspace(log10(smallest),log10(largest),50);
    profileTheta = sort(unique([baseTheta,opt.phi-offsets,opt.phi+offsets]));
else
    profileTheta = baseTheta;
end
profileTheta = profileTheta(profileTheta>0 & profileTheta<pi & ...
    abs(profileTheta-opt.phi)>1e-13);
profileTheta = profileTheta(:);
p = c+a*[sin(profileTheta),zeros(numel(profileTheta),1), ...
    cos(profileTheta)];
profile = zeros(numel(profileTheta),1);
for first = 1:opt.batchSize:numel(profileTheta)
    rows = first:min(first+opt.batchSize-1,numel(profileTheta));
    u = stokes_mfs_velocity_matrix( ...
        exteriorBasis,p(rows,:))*(solution.lambda/opt.muOut);
    u = reshape(u,3,[]).'-solution.U(:).';
    tangent = [cos(profileTheta(rows)),zeros(numel(rows),1), ...
        -sin(profileTheta(rows))];
    profile(rows) = sum(u.*tangent,2);
end
mobile = profileTheta(:)>opt.phi;
s = a*(profileTheta(mobile)-opt.phi);
mobileVelocity = profile(mobile);
diagnostic.fitC = NaN;
diagnostic.fitWindow = [NaN NaN];
diagnostic.profileTheta = profileTheta(:);
diagnostic.profileVelocity = profile;
diagnostic.closestMobileS = min(s);
diagnostic.closestMobileVelocity = mobileVelocity(1);
if ~isempty(solution.grid.ringOffsets)
    lower = 4*a*solution.grid.ringOffsets(end)/4;
    upper = a*opt.edgeBandAngle/2;
    fit = s>=lower & s<=upper;
    if nnz(fit)>=2
        diagnostic.fitC = sqrt(s(fit))\mobileVelocity(fit);
        diagnostic.fitWindow = [lower upper];
    end
end
profileFigure = figure('Name','Stagnant-cap meridional profile', ...
    'Color','w','Visible',visibility);
tiledlayout(profileFigure,1,3);
nexttile;
plot(profileTheta,profile,'b-','LineWidth',1.4); hold on;
xline(opt.phi,'k--');
plot(opt.phi,0,'ko','MarkerFaceColor','k');
xlabel('\theta'); ylabel('u_\theta-U\cdot e_\theta');
title('Surface profile'); grid on;
nexttile;
numericalHandle = plot(s,mobileVelocity,'bo-','MarkerSize',3); hold on;
if isfinite(diagnostic.fitC)
    referenceS = linspace(0,diagnostic.fitWindow(2),80);
    squareRootHandle = plot(referenceS,diagnostic.fitC*sqrt(referenceS), ...
        'r--','LineWidth',1.4);
    xline(diagnostic.fitWindow(1),'k:');
end
zeroHandle = plot(0,0,'ko','MarkerFaceColor','k');
xlim([0,max(a*opt.edgeBandAngle,eps)]);
if isfinite(diagnostic.fitC)
    legend([numericalHandle,squareRootHandle,zeroHandle], ...
        {'computed','C\surd s fit','exact edge limit'},'Location','best');
end
xlabel('s=a(\theta-\phi)'); ylabel('relative u_\theta');
title('Mobile-side edge'); grid on;
nexttile;
plot(s,mobileVelocity./sqrt(s),'bo-','MarkerSize',3);
if isfinite(diagnostic.fitC), yline(diagnostic.fitC,'r--'); end
xlim([0,max(a*opt.edgeBandAngle,eps)]);
xlabel('s=a(\theta-\phi)'); ylabel('u_\theta/\surd s');
title('Square-root ratio'); grid on;
sgtitle('Exact law: u_\theta\sim C\surd s, derivative \sim C/(2\surd s); finite MFS trace shown');
end

function self_test()
%SELF_TEST Check the plotted meridional trace against a direct kernel value.
o = stagnant_cap_options(struct('sourceRequest',50,'ringLevels',2, ...
    'plotThetaCount',24,'plotAzimuthCount',24));
g = stagnant_cap_geometry(o);
s.grid = g;
s.lambda = zeros(3*g.N,1);
s.lambda(3) = 1;
s.U = zeros(3,1);
[f1,f2,d] = plot_stagnant_cap_velocity(s,o,'off');
assert(isgraphics(f1) && isgraphics(f2));
theta = d.profileTheta(1);
point = o.centre+o.radius*[sin(theta),0,cos(theta)];
u = drop_velocity_block(g.innerSources,point,o.muOut)*s.lambda;
reference = u.'*[cos(theta);0;-sin(theta)];
assert(abs(reference-d.profileVelocity(1))<1e-12);
close(f1); close(f2);
fprintf('plot_stagnant_cap_velocity self-test passed.\n');
end
