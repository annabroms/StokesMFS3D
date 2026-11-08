function result = ring_family_case(Y0,id0,Yc,idc,centre,a,mu,Fg,families, ...
    geom,op,val,audit,Uexact,beta,tol,mode,phi)
%RING_FAMILY_CASE Solve and audit one mixed-family axisymmetric ring case.
% Inputs Y0/Yc: N0/Nc-by-3 sources; id0/idc: ring membership; centre: 1-by-3;
% a,mu: radius and viscosity; Fg,Uexact: 3-by-1 axial vectors; families:
% ordered ring_family_list tokens. geom and op are shared surface geometry
% and frame maps from stagnant_cap_geometry and bubble_boundary_operators.
% val and audit contain independent Mv-by-3 surface nodes, unit normals and
% cap/mobile masks. beta is velocity/force, tol is the relative SVD cutoff.
% Coefficients are block-major, then ring-major [rho;psi;z]; C reconstructs
% node-interleaved Cartesian vectors. Boundary velocities are normalized by
% |Uexact| and tangential traction by mu*|Uexact|/a. Full vector velocity
% holds on the cap; mobile rows enforce normal velocity and tangent traction.
% result contains audit maps, max residuals, force/torque, swirl and per-family
% maximum per-node coefficient divided by a documented characteristic scale.
% Optional mode is "full" (default) or "swirl_free"; see ring_mixed_basis.
% Optional phi is the physical cap angle in radians used for regional
% diagnostics; by default it is inferred from the first source as before.
if nargin==0
    self_test();
    return
end
if nargin<17, mode = "full"; end
if isempty(families)
    familyList = strings(1,0);
    result.name = "GL only";
else
    familyList = ring_family_list(families);
    result.name = "GL + "+join(familyList,"+");
end
[basis,C,w,k,moments,selected] = ring_mixed_basis(Y0,id0,Yc,idc,centre,a,familyList,mode);
Mc = nnz(geom.cap); Mm = nnz(geom.mobile);
G = stokes_mfs_velocity_matrix(basis,geom.points)/mu;
T = stokes_mfs_traction_matrix(basis,geom.points(geom.mobile,:),geom.normals(geom.mobile,:));
capNodes = find(geom.cap);
cap = reshape([3*capNodes-2,3*capNodes-1,3*capNodes].',[],1);
mobile = op.mobileVectorRows;
A = [G(cap,:)*C;op.Nm*G(mobile,:)*C;op.Dm*T*C];
V = [op.KMc(:,3);op.Nm*op.KMm(:,3);zeros(2*Mm,1)];
Uref = norm(Uexact); tractionScale = mu*Uref/a;
rowScale = [ones(3*Mc+Mm,1)/Uref;ones(2*Mm,1)/tractionScale];
familyNames = ["GL S",familyList];
familyScale = zeros(1,numel(familyNames));
familyScale(1) = norm(Fg)/size(Y0,1);
colScale = repmat(familyScale(1),3*max(id0),1);
for j = 1:numel(familyList)
    switch familyList(j)
        case {"R","T","T_t","T_n","T_n2"}
            dimensionalScale = a;
        case "D"
            dimensionalScale = a^2;
        otherwise
            dimensionalScale = 1;
    end
    familyScale(j+1) = norm(Fg)*dimensionalScale/size(Yc,1);
    colScale = [colScale;repmat(familyScale(j+1),3*max(idc),1)]; %#ok<AGROW>
end
colScale = colScale(selected);
fprintf('\n%s: A %d x %d; per-family coefficient scales:',result.name,size(A));
fprintf(' %.3g',familyScale); fprintf('\n');
[~,b,Uz,info] = ring_complete_solve(A,V,w,k,k*Fg(3),1,beta,rowScale,colScale,tol);
lambda = C*b;
solution = struct('basis',basis,'lambda',lambda,'U',[0;0;Uz]);
evalOpt = struct('centre',centre,'radius',a,'muOut',mu,'batchSize',128);
if nargin<18
    phi = acos((Yc(1,3)-centre(3))/norm(Yc(1,:)-centre));
end
[result.residual,~,~,~,~] = check_grid(solution,val,evalOpt,Uref,tractionScale,centre,phi);
[result.auditResidual,result.mapPoints,result.mapResiduals,result.swirl,result.regional] = ...
    check_grid(solution,audit,evalOpt,Uref,tractionScale,centre,phi);
moment = moments*b;
result.forceResidual = norm(moment(1:3)-Fg)/norm(Fg);
result.torqueResidual = norm(moment(4:6))/(a*norm(Fg));
result.Uerror = norm(solution.U-Uexact)/Uref;
result.U = solution.U; result.rank = info.rank; result.unknowns = size(A,2);
result.familyNames = familyNames; result.coefficientScale = familyScale;
result.coefficientMax = zeros(size(familyScale));
offset = 0;
for j = 1:numel(basis)
    count = 3*size(basis(j).points,1);
    vectors = reshape(lambda(offset+(1:count)),3,[]);
    result.coefficientMax(j) = max(vecnorm(vectors,2,1))/familyScale(j);
    offset = offset+count;
end
end

function [maxima,points,residuals,swirl,regional] = check_grid(solution,grid,opt,Uref,tractionScale,centre,phi)
%CHECK_GRID Evaluate the physical fields and region-specific residual maps.
field = bubble_evaluate(solution,grid.points,opt,grid.normals);
cap = grid.cap; mobile = grid.mobile;
points = {grid.points(cap,:),grid.points(mobile,:),grid.points(mobile,:)};
residuals = {vecnorm(field.relativeVelocity(cap,:),2,2)/Uref, ...
    abs(field.normalVelocity(mobile))/Uref, ...
    vecnorm(field.tangentTraction(mobile,:),2,2)/tractionScale};
maxima = nan(1,3);
for j = 1:3
    if ~isempty(residuals{j}), maxima(j) = max(residuals{j}); end
end
psi = atan2(grid.points(:,2)-centre(2),grid.points(:,1)-centre(1));
ephi = [-sin(psi),cos(psi),zeros(numel(psi),1)];
swirl = max(abs(sum(field.velocity.*ephi,2)))/Uref;
theta = acos(max(-1,min(1,grid.normals(:,3))));
regional = nan(3,3);
regions = {abs(theta-phi)<=0.1,abs(theta-phi)>0.1,theta<0.1 | theta>pi-0.1};
conditionMasks = {cap,mobile,mobile};
for j = 1:3
    for q = 1:3
        selected = regions{q}(conditionMasks{j});
        if any(selected), regional(j,q) = max(residuals{j}(selected)); end
    end
end
end

function self_test()
%SELF_TEST Exercise a small mixed solve and verify force and finite fields.
n = 8; psi = (0:n-1)'*2*pi/n; centre = [0 0 0];
Y0 = [0.35*cos(psi),0.35*sin(psi),0.1*ones(n,1)];
Yc = [0.5*cos(psi),0.5*sin(psi),0.5*ones(n,1)];
dirs = [eye(3);-eye(3);1 1 1;1 -1 1;-1 1 1;-1 -1 1];
dirs = dirs./vecnorm(dirs,2,2); M = size(dirs,1);
cap = dirs(:,3)>0.5;
geom = struct('sourcePoints',Y0,'baseSourceCount',n,'coefficientCount',3*n, ...
    'collocationCount',M,'points',dirs,'normals',dirs,'cap',cap,'mobile',~cap);
op = bubble_boundary_operators(geom,struct('n_level',0,'n_az_source',0,'muOut',1));
grid = struct('points',dirs,'normals',dirs,'cap',cap,'mobile',~cap);
Fg = [0;0;-1]; Uexact = Fg/(5*pi);
sets = {strings(1,0),"S","SD","SR","ST_t","ST_tT_n","SQ_s","ST_n2","SS_y","ST_nS_y"};
for j = 1:numel(sets)
    r = ring_family_case(Y0,ones(n,1),Yc,ones(n,1),centre,1,1,Fg, ...
        sets{j},geom,op,grid,grid,Uexact,1/(6*pi),1e-10);
    assert(all(isfinite([r.U;r.residual(:);r.auditResidual(:)])));
    assert(r.forceResidual<1e-12 && r.rank>0);
    assert(all(r.auditResidual==r.residual));
    r = ring_family_case(Y0,ones(n,1),Yc,ones(n,1),centre,1,1,Fg, ...
        sets{j},geom,op,grid,grid,Uexact,1/(6*pi),1e-10,"swirl_free");
    assert(all(isfinite([r.U;r.residual(:);r.auditResidual(:)])));
    assert(r.forceResidual<1e-12 && r.rank>0);
    if isempty(sets{j})
        expected = 2;
    else
        families = ring_family_list(sets{j});
        expected = 2+2*numel(families)-nnz(ismember(families,["R","T_n2"]));
    end
    assert(r.unknowns==expected);
end
% The GL-only three-component basis must preserve the no-swirl solution
% on this reflection-symmetric surface grid and axial forcing.
angle = 0.19;
Q = [cos(angle),-sin(angle),0;sin(angle),cos(angle),0;0,0,1];
check = grid;
check.points = dirs*Q.';
check.normals = check.points;
r = ring_family_case(Y0,ones(n,1),Yc,ones(n,1),centre,1,1,Fg, ...
    strings(1,0),geom,op,check,check,Uexact,1/(6*pi),1e-10);
[C2,W2,K2] = ring_cylindrical_maps(ones(n,1),psi);
mc = find(cap);
capRows = reshape([3*mc-2,3*mc-1,3*mc].',[],1);
G2 = op.G*C2;
A2 = [G2(capRows,:);op.Nm*G2(op.mobileVectorRows,:);op.Dm*op.Tm*C2];
V2 = [op.KMc(:,3);op.Nm*op.KMm(:,3);zeros(2*nnz(~cap),1)];
scale = norm(Uexact);
[~,b2,U2] = ring_complete_solve(A2,V2,W2,K2,K2*Fg(3)/n,n, ...
    1/(6*pi),ones(3*M,1)/scale,1/n,1e-10);
assert(abs(U2-r.U(3))<1e-11);
r2 = ring_family_case(Y0,ones(n,1),Yc,ones(n,1),centre,1,1,Fg, ...
    strings(1,0),geom,op,check,check,Uexact,1/(6*pi),1e-10,"swirl_free");
assert(abs(U2-r2.U(3))<1e-11 && r2.unknowns==2);
inferredPhi = acos(Yc(1,3)/norm(Yc(1,:)));
r3 = ring_family_case(Y0,ones(n,1),Yc,ones(n,1),centre,1,1,Fg, ...
    strings(1,0),geom,op,check,check,Uexact,1/(6*pi),1e-10,"swirl_free",inferredPhi);
assert(isequaln(r2.regional,r3.regional) && norm(r2.U-r3.U)<1e-14);
reference = struct('basis',op.basis,'lambda',C2*b2,'U',[0;0;U2]);
evalOpt = struct('centre',centre,'radius',1,'muOut',1,'batchSize',4);
f = bubble_evaluate(reference,check.points,evalOpt,check.normals);
referenceResiduals = {vecnorm(f.relativeVelocity(cap,:),2,2)/scale, ...
    abs(f.normalVelocity(~cap))/scale, ...
    vecnorm(f.tangentTraction(~cap,:),2,2)/scale};
for j = 1:3
    assert(norm(referenceResiduals{j}-r.mapResiduals{j})<1e-11);
    assert(norm(referenceResiduals{j}-r2.mapResiduals{j})<1e-11);
end
fprintf('ring_family_case self-test passed.\n');
end
