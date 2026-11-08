%% Ring-constant Cartesian Stokeslets on a sphere (note Section 10, rev. 8)
% One sphere in free space, prescribed force Fg along z, unknown translation U.
% Stokeslets only on a concentric interior Gauss-Legendre proxy sphere; every
% node of ring r carries the same Cartesian vector a_r (no azimuthal rotation).
% Exactly two solves on identical geometry: phi = 0 (clean bubble) and
% phi = pi (rigid sphere). Run after setpath. Nodal vectors are interleaved
% [x1;y1;z1;x2;...]; ring vectors a are ring-major [a_1;...;a_R].

%% Parameters
a      = 1;            % physical sphere radius
mu     = 1;            % exterior viscosity
centre = [0 0 0];
Fg     = [0;0;-1];     % prescribed force on the body
ap     = 0.6*a;        % proxy sphere radius (sources only here)
Nu     = 32;           % azimuthal resolution scale of the GL proxy sphere
Nv     = 16;           % number of GL latitude rings R (rounded up to even)
Mreq   = 600;          % requested collocation count (spherical design)
validationMultiplier = 5;   % validation grid has >= this many times M nodes
svdTol = 1e-10;        % relative truncation threshold of the scaled matrix
quiverFraction = 0.3;  % longest arrow = this fraction of the proxy radius
ringBasis = 'cylindrical';  % 'cartesian' (Section 10) or 'cylindrical' (a_rho,a_z per ring)
beta   = 1/(6*pi*mu*a);     % velocity/force conversion scale, not the mobility
Uref   = [norm(Fg)/(4*pi*mu*a), norm(Fg)/(6*pi*mu*a)];  % bubble, rigid
Uexact = [Fg/(4*pi*mu*a), Fg/(6*pi*mu*a)];
caseNames = ["clean bubble (phi=0)","rigid sphere (phi=pi)"];

%% Geometry
[Y,ringId,nr,theta] = ring_gl_sphere(ap,centre,Nu,Nv);
R = numel(nr); N = size(Y,1);
cyl = strcmp(ringBasis,'cylindrical');
if cyl      % unknowns [a_rho,r;a_z,r]; only the z force and z translation remain
    [C,W,KR] = ring_cylindrical_maps(ringId,atan2(Y(:,2)-centre(2),Y(:,1)-centre(1)));
    nu = 2; Fc = Fg(3); zc = 3;
else
    [C,W,KR] = ring_constant_maps(ringId);
    nu = 3; Fc = Fg; zc = 1:3;
end
dirs = stagnant_cap_design(Mreq);               % collocation: spherical design
M = size(dirs,1);
xc = centre+a*dirs; nc = dirs;
gridOpt = struct('radius',a,'centre',centre,'muOut',mu,'phi',0, ...
    'sourceRequest',50,'validationMultiplier',validationMultiplier);
val = stagnant_cap_validation_grid(gridOpt,M,1);   % independent rotated grid
fprintf('Rings R=%d, sources N=%d, collocation M=%d, validation Mv=%d, scalar unknowns %d (%s).\n', ...
    R,N,M,val.count,nu*R,ringBasis);

%% Assembly
% Nodal operators from the existing routine (all targets mobile, nodal basis).
geom = struct('sourcePoints',Y,'baseSourceCount',N,'coefficientCount',3*N, ...
    'collocationCount',M,'points',xc,'normals',nc, ...
    'cap',false(M,1),'mobile',true(M,1));
op = bubble_boundary_operators(geom,struct('n_level',0,'n_az_source',0,'muOut',mu));
GR = op.G*C;  TR = op.Tm*C;     % sum nodal columns ring by ring, eq. (65)
KM = op.KM;
Nm = op.Nm;  Dm = op.Dm;        % normal (M x 3M) and tangent (2M x 3M) maps
fprintf('Nodal G %dx%d -> ring GR %dx%d; TR %dx%d.\n', ...
    size(op.G),size(GR),size(TR));
a0 = KR*Fc/N;                    % completion: every Stokeslet carries Fg/N
% Case 1 rows: normal velocity (M), tangential traction (2M). Case 2: velocity (3M).
Aset = {[Nm*GR;Dm*TR],GR};
Vset = {[Nm*KM(:,zc);zeros(2*M,numel(zc))],KM(:,zc)};
tractionScale = mu*Uref(1)/a;
rowScale = {[ones(M,1)/Uref(1);ones(2*M,1)/tractionScale],ones(3*M,1)/Uref(2)};

%% Solution
sol = cell(1,2);
for k = 1:2
    fprintf('\n%s: A %dx%d, V %dx%d\n',caseNames(k),size(Aset{k}),size(Vset{k}));
    [alpha,ar,U,info] = ring_complete_solve(Aset{k},Vset{k},W,KR,a0,N,beta, ...
        rowScale{k},norm(Fg),svdTol);
    Uf = zeros(3,1); Uf(zc) = U;        % full translation vector
    sol{k} = struct('alpha',alpha,'a',ar,'U',Uf,'info',info);
end

%% Validation on the independent grid
evalOpt = struct('centre',centre,'radius',a,'muOut',mu,'batchSize',256);
res = cell(1,3);                 % bubble normal, bubble traction, rigid velocity
fprintf('\n%-24s %-12s %-12s %-12s %-12s\n','case','U_z','rel. U err','force res','torque/(|Fg|a)');
for k = 1:2
    lambda = C*sol{k}.a;
    s = struct('basis',op.basis,'lambda',lambda,'U',sol{k}.U);
    f = bubble_evaluate(s,val.points,evalOpt,val.normals);
    forceRes = norm(sum(reshape(lambda,3,[]),2)-Fg)/norm(Fg);
    lam = reshape(lambda,3,[]).';
    torque = -sum(cross(Y-centre,lam,2),1);
    if k == 1
        res{1} = abs(f.normalVelocity)/Uref(1);
        res{2} = vecnorm(f.tangentTraction,2,2)/tractionScale;
        maxima = [max(res{1}),max(res{2})];
    else
        res{3} = vecnorm(f.relativeVelocity,2,2)/Uref(2);
        maxima = max(res{3});
    end
    fprintf('%-24s %-12.6f %-12.3e %-12.3e %-12.3e\n',caseNames(k), ...
        sol{k}.U(3),norm(sol{k}.U-Uexact(:,k))/norm(Uexact(:,k)), ...
        forceRes,norm(torque)/(norm(Fg)*a));
    fprintf('   exact U_z %.6f; transverse |U_xy| %.2e; max normalized residuals:',Uexact(3,k), ...
        norm(sol{k}.U(1:2)));
    fprintf(' %.3e',maxima); fprintf('\n');
end
fprintf('Bubble residuals: normal velocity / (|Fg|/(4 pi mu a)), traction / (mu U/a).\n');
fprintf('Rigid residual: |u-U| / (|Fg|/(6 pi mu a)).\n');

%% Plotting
ring_plot_residuals(val.points,res, ...
    {'Bubble |(u-U).n| / U_{bubble}','Bubble tangential traction / (\mu U_{bubble}/a)', ...
     'Rigid |u-U| / U_{rigid}'});
first = [1;cumsum(nr(1:end-1))+1];            % psi = 0 node of each ring
lamRep = cellfun(@(s) reshape(C*s.a,3,[]).',sol,'UniformOutput',false);   % per-node vectors
arrowScale = quiverFraction*ap/max(cellfun(@(l) max(vecnorm(l(first,:),2,2)),lamRep));
figure('Color','w','Position',[60 100 1000 480]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for k = 1:2
    ring_plot_sources(nexttile(tl),Y(first,:),lamRep{k}(first,:), ...
        centre,a,ap,arrowScale,caseNames(k));
end
