%% Ring-constant sphere with a stagnant cap: sweep in the cap angle phi
% Same geometry, ring maps and force completion as oct29_ring_constant_sphere.m,
% but the cap theta < phi (polar angle from +z, about the sphere centre) is
% stagnant for every phi in the sweep. One solve per phi, on identical sources.
%   cap     (theta < phi): full velocity u = U               (3 rows per node)
%   mobile  (theta >= phi): normal velocity u.n = U.n and zero tangential
%                           traction                          (1 + 2 rows per node)
% phi = 0 is the clean bubble and phi = pi the rigid sphere. Exact translation:
% U = Fg/zeta(phi), zeta from stagnant_cap_exact_drag (Sadhal-Johnson).
% Run after setpath. Nodal vectors are interleaved [x1;y1;z1;x2;...]; ring
% unknowns are ring-major, see ring_constant_maps / ring_cylindrical_maps.

%% Parameters
phis   = [0 1/12 1/6 1/3 1/2 2/3 5/6 11/12 1]*pi;   % cap angles to sweep
a      = 1;            % physical sphere radius
mu     = 1;            % exterior viscosity
centre = [0 0 0];
Fg     = [0;0;-1];     % prescribed force on the body
ap     = 0.6*a;        % proxy sphere radius (sources only here)
Nu     = 40;           % azimuthal resolution scale of the GL proxy sphere
Nv     = 30;           % number of GL latitude rings R (rounded up to even)
Mreq   = 600;          % requested base collocation count (spherical design)
ringLevels = 8;        % extra collocation rings on each side of the cap edge
clusterSigma = 3;      % tapered clustering strength, offsets e^{-sigma(sqrt(L)-sqrt(L-k))}
phiShow = [1/6 1/2 5/6]*pi;  % cap angles with surface residual maps (members of phis)
validationMultiplier = 5;   % validation grid has >= this many times M nodes
edgeOffset = 1e-5;     % finest validation ring offset from the cap edge (rad)
svdTol = 1e-10;        % relative truncation threshold of the scaled matrix
ringBasis = 'cylindrical'; % 'cartesian' (Section 10.1) or 'cylindrical' (a_rho,a_z)
beta   = 1/(6*pi*mu*a);    % velocity/force conversion scale, not the mobility

%% Geometry (shared by all phi)
[Y,ringId,nr] = ring_gl_sphere(ap,centre,Nu,Nv);
R = numel(nr); N = size(Y,1);
cyl = strcmp(ringBasis,'cylindrical');
if cyl      % unknowns [a_rho,r;a_z,r]; only the z force and z translation remain
    [C,W,KR] = ring_cylindrical_maps(ringId,atan2(Y(:,2)-centre(2),Y(:,1)-centre(1)));
    nu = 2; Fc = Fg(3); zc = 3;
else
    [C,W,KR] = ring_constant_maps(ringId);
    nu = 3; Fc = Fg; zc = 1:3;
end
a0 = KR*Fc/N;                                   % completion: every Stokeslet carries Fg/N
% Collocation per phi: spherical design plus tapered edge rings, reusing
% stagnant_cap_geometry (its proxy sources are not used here).
colloc = struct('radius',a,'centre',centre,'muOut',mu,'sourceRequest',50, ...
    'collocationRatio',Mreq/size(stagnant_cap_design(50),1),'ringLevels',ringLevels, ...
    'collocClustering',"tapered",'collocClusterSigma',clusterSigma);
fprintf('Rings R=%d, sources N=%d, scalar unknowns %d (%s).\n',R,N,nu*R,ringBasis);

% Nodal kernels do not depend on phi except through which rows are cap rows,
% so the operators are rebuilt per phi from the same sources and nodes.
nPhi = numel(phis);
out = struct('Uz',zeros(nPhi,1),'UzExact',zeros(nPhi,1),'Uerr',zeros(nPhi,1), ...
    'capVel',nan(nPhi,1),'mobNormal',nan(nPhi,1),'mobTraction',nan(nPhi,1), ...
    'force',zeros(nPhi,1),'rank',zeros(nPhi,1),'Mv',zeros(nPhi,1), ...
    'maxStrength',zeros(nPhi,1),'M',zeros(nPhi,1));
maps = cell(nPhi,1);             % per-phi validation points/residuals for the maps
evalOpt = struct('centre',centre,'radius',a,'muOut',mu,'batchSize',256);

%% Sweep: assembly, solution, validation
for k = 1:nPhi
    phi = phis(k);
    cg = stagnant_cap_geometry(setfield(colloc,'phi',phi)); %#ok<SFLD>
    M = cg.collocationCount; cap = cg.cap;
    geom = struct('sourcePoints',Y,'baseSourceCount',N,'coefficientCount',3*N, ...
        'collocationCount',M,'points',cg.points,'normals',cg.normals, ...
        'cap',cap,'mobile',cg.mobile);
    op = bubble_boundary_operators(geom,struct('n_level',0,'n_az_source',0,'muOut',mu));
    Gm = op.G(op.mobileVectorRows,:);           % mobile-row velocity (3Mm x 3N)
    % Rows: cap velocity (3Mc), mobile normal velocity (Mm), mobile traction (2Mm).
    A = [op.Gc*C; op.Nm*Gm*C; op.Dm*op.Tm*C];
    V = [op.KMc; op.Nm*op.KMm; zeros(2*nnz(~cap),3)];
    V = V(:,zc);
    zeta = stagnant_cap_exact_drag(phi,0,mu,a);  % bubble: inner viscosity zero
    Uexact = Fg/zeta; Uref = norm(Uexact); tractionScale = mu*Uref/a;
    rowScale = [ones(3*nnz(cap)+nnz(~cap),1)/Uref; ones(2*nnz(~cap),1)/tractionScale];
    fprintf('\nphi = %.4f pi: collocation %d (cap %d, mobile %d), A %dx%d\n', ...
        phi/pi,M,nnz(cap),nnz(~cap),size(A));
    [~,ar,U,info] = ring_complete_solve(A,V,W,KR,a0,N,beta,rowScale,norm(Fg),svdTol);
    Uf = zeros(3,1); Uf(zc) = U;

    % Independent validation grid with edge rings for this phi.
    gridOpt = struct('radius',a,'centre',centre,'muOut',mu,'phi',phi, ...
        'sourceRequest',50,'validationMultiplier',validationMultiplier);
    if phi>0 && phi<pi
        val = stagnant_cap_validation_grid(gridOpt,M,edgeOffset);
    else
        val = stagnant_cap_validation_grid(gridOpt,M,1);
    end
    lambda = C*ar;
    f = bubble_evaluate(struct('basis',op.basis,'lambda',lambda,'U',Uf), ...
        val.points,evalOpt,val.normals);
    vc = val.cap; vm = val.mobile;
    rCap = vecnorm(f.relativeVelocity(vc,:),2,2)/Uref;
    rNormal = abs(f.normalVelocity(vm))/Uref;
    rTraction = vecnorm(f.tangentTraction(vm,:),2,2)/tractionScale;
    out.capVel(k) = max([rCap;NaN]); out.mobNormal(k) = max([rNormal;NaN]);
    out.mobTraction(k) = max([rTraction;NaN]);
    maps{k} = struct('points',{{val.points(vc,:),val.points(vm,:),val.points(vm,:)}}, ...
        'res',{{rCap,rNormal,rTraction}});
    out.maxStrength(k) = max(vecnorm(reshape(lambda,3,[]),2,1));   % per-node |lambda|
    out.M(k) = M;
    out.Uz(k) = Uf(3); out.UzExact(k) = Uexact(3);
    out.Uerr(k) = norm(Uf-Uexact)/norm(Uexact);
    out.force(k) = norm(sum(reshape(lambda,3,[]),2)-Fg)/norm(Fg);
    out.rank(k) = info.rank; out.Mv(k) = val.count;
end

%% Report
fprintf('\n%-9s %-10s %-10s %-10s %-11s %-11s %-11s %-9s\n','phi/pi','U_z','exact U_z', ...
    'rel U err','cap |u-U|','mob u.n','mob traction','force res');
for k = 1:nPhi
    fprintf('%-9.4f %-10.6f %-10.6f %-10.2e %-11.2e %-11.2e %-11.2e %-9.1e\n', ...
        phis(k)/pi,out.Uz(k),out.UzExact(k),out.Uerr(k),out.capVel(k), ...
        out.mobNormal(k),out.mobTraction(k),out.force(k));
end
fprintf(['Normalizations: velocities by |U_exact(phi)|, traction by mu|U_exact(phi)|/a; ' ...
    'NaN = region empty at this phi.\n']);

%% Plotting
figure('Color','w','Position',[60 100 1500 520]);
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl);
semilogy(ax,phis/pi,out.capVel,'o-',phis/pi,out.mobNormal,'s-', ...
    phis/pi,out.mobTraction,'^-','LineWidth',1.2);
grid(ax,'on'); xlabel(ax,'\phi/\pi');
ylabel(ax,'max normalized residual (validation grid)');
legend(ax,{'cap |u-U|/U_{ex}','mobile |(u-U).n|/U_{ex}', ...
    'mobile tangential traction/(\mu U_{ex}/a)'},'Location','southoutside','FontSize',8);
title(ax,sprintf('Boundary residuals, %s ring basis, R=%d',ringBasis,R));
ax = nexttile(tl);
semilogy(ax,phis/pi,max(out.Uerr,eps),'o-','LineWidth',1.2);
grid(ax,'on'); xlabel(ax,'\phi/\pi');
ylabel(ax,'|U-U_{ex}|/|U_{ex}|');
title(ax,'Translation error vs Sadhal-Johnson U = F_g/\zeta(\phi)');
ax = nexttile(tl);
semilogy(ax,phis/pi,out.maxStrength*N/norm(Fg),'o-','LineWidth',1.2);
grid(ax,'on'); xlabel(ax,'\phi/\pi');
ylabel(ax,'max_{rj} |\lambda_{rj}| / (|F_g|/N)');
title(ax,'Largest per-node Stokeslet strength');

% Surface maps of where each condition's residual is largest (check nodes).
for phi = phiShow
    m = maps{abs(phis-phi)<1e-12};
    ring_plot_residuals(m.points,m.res, ...
        {sprintf('phi=%.3g pi cap |u-U|/U_ex',phi/pi), ...
         sprintf('phi=%.3g pi mobile |(u-U).n|/U_ex',phi/pi), ...
         sprintf('phi=%.3g pi mobile tang. traction/(mu U_ex/a)',phi/pi)});
end
