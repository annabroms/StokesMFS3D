function system = stagnant_cap_system(geometry,opt,drive)
%STAGNANT_CAP_SYSTEM Assemble four boundary row blocks for mixed source bases.
% Rows are [continuity(3M); normal(M); cap tangent(2Mc);
% mobile traction(2Mm)]. Columns are [exterior coefficients;interior
% coefficients], with interleaved x,y,z components in each source block.
% All tractions use n toward exterior. Projection preserves the full
% exterior source-force map, including every added Stokeslet.
if nargin == 0
    self_test();
    return
end
opt = stagnant_cap_options(opt);
M = geometry.collocationCount;
if isfield(geometry,'exteriorBasis')
    exteriorBasis = geometry.exteriorBasis;
    interiorBasis = geometry.interiorBasis;
else
    block = struct('family',"S",'points',geometry.innerSources, ...
        'normalDirections',zeros(0,3),'scale',1);
    exteriorBasis = block;
    block.points = geometry.outerSources;
    interiorBasis = block;
end
outerColumns = 3*sum(arrayfun( ...
    @(b) size(b.points,1),exteriorBasis));
innerColumns = 3*sum(arrayfun( ...
    @(b) size(b.points,1),interiorBasis));
cap = find(geometry.cap);
mobile = find(geometry.mobile);
Mc = numel(cap);
Mm = numel(mobile);
assert(Mc+Mm==M && isequal(size(geometry.normals),[M 3]));
Nmap = zeros(M,3*M);
Dmap = zeros(2*M,3*M);
reference = repmat([1 0 0],M,1);
reference(abs(geometry.normals(:,1))>0.9,:) = ...
    repmat([0 1 0],nnz(abs(geometry.normals(:,1))>0.9),1);
t1 = cross(geometry.normals,reference,2);
t1 = t1./vecnorm(t1,2,2);
t2 = cross(geometry.normals,t1,2);
for p = 1:M
    columns = 3*p-2:3*p;
    Nmap(p,columns) = geometry.normals(p,:);
    Dmap(2*p-1:2*p,columns) = [t1(p,:);t2(p,:)];
end
capTangentRows = reshape([2*cap-1,2*cap].',[],1);
mobileTangentRows = reshape([2*mobile-1,2*mobile].',[],1);
Dc = Dmap(capTangentRows,:);
Dm = Dmap(mobileTangentRows,:);
KM = kron(ones(M,1),eye(3));
G = stokes_mfs_velocity_matrix(exteriorBasis,geometry.points)/opt.muOut;
E = stokes_mfs_velocity_matrix(interiorBasis,geometry.points)/opt.muIn;
T = stokes_mfs_traction_matrix( ...
    exteriorBasis,geometry.points,geometry.normals);
Q = stokes_mfs_traction_matrix( ...
    interiorBasis,geometry.points,geometry.normals);
assert(isequal(size(G),[3*M outerColumns]) && ...
    isequal(size(E),[3*M innerColumns]) && ...
    isequal(size(T),[3*M outerColumns]) && ...
    isequal(size(Q),[3*M innerColumns]));
body.basis = exteriorBasis;
[moments,~] = image_mobility_maps(body,opt.centre);
forceMap = moments(1:3,:);
forceCompletion = forceMap'/(forceMap*forceMap');
assert(norm(forceMap*forceCompletion-eye(3),'fro')<1e-11);
zeroNormal = zeros(M,innerColumns);
zeroCap = zeros(2*Mc,innerColumns);
kind = string(drive.kind);
if kind=="force"
    Fg = drive.Fg(:);
    assert(numel(Fg)==3 && isreal(Fg) && all(isfinite(Fg)) && ...
        norm(Fg)>0 && norm(Fg(1:2))<1e-13*norm(Fg), ...
        'Force must be nonzero and axial.');
    lambda0 = forceCompletion*Fg;
    % Apply Pi_F=I-C_F*H_F without forming the full projector.
    Gp = G-(G*forceCompletion)*forceMap;
    Tp = T-(T*forceCompletion)*forceMap;
    beta = opt.completionScale;
    system.matrix = [Gp,-E; ...
        Nmap*Gp+beta*(Nmap*KM)*forceMap,zeroNormal; ...
        Dc*Gp+beta*(Dc*KM)*forceMap,zeroCap; ...
        Dm*Tp,-Dm*Q];
    system.rhs = -[G*lambda0;Nmap*G*lambda0; ...
        Dc*G*lambda0;Dm*T*lambda0];
    system.velocityScale = norm(Fg)/(6*pi*opt.muOut*opt.radius);
    system.lambda0 = lambda0;
    system.force = Fg;
elseif kind=="speed"
    Uz = drive.Uz;
    assert(isscalar(Uz) && isreal(Uz) && isfinite(Uz) && Uz~=0);
    U = [0;0;Uz];
    system.matrix = [G,-E; ...
        Nmap*G,zeroNormal; ...
        Dc*G,zeroCap; ...
        Dm*T,-Dm*Q];
    system.rhs = [zeros(3*M,1);Nmap*KM*U; ...
        Dc*KM*U;zeros(2*Mm,1)];
    system.velocityScale = abs(Uz);
    system.lambda0 = zeros(outerColumns,1);
    system.force = [];
else
    error('stagnant_cap_system:BadMode','Use force or speed mode.');
end
system.kind = kind;
system.KN = forceMap.'; % Generalized source-force transpose.
system.forceMap = forceMap;
system.momentMap = moments;
system.forceCompletion = forceCompletion;
system.exteriorBasis = exteriorBasis;
system.interiorBasis = interiorBasis;
system.outerColumns = outerColumns;
system.innerColumns = innerColumns;
system.M = M;
system.Mc = Mc;
system.Mm = Mm;
system.tangentOne = t1;
system.tangentTwo = t2;
system.tractionScale = opt.muOut*system.velocityScale/opt.radius;
assert(isequal(size(system.matrix),[6*M,outerColumns+innerColumns]));
assert(numel(system.rhs)==6*M);
end

function self_test()
%SELF_TEST Check the old all-mobile block and both translation completions.
o = stagnant_cap_options(struct('sourceRequest',50,'phi',0, ...
    'muOut',1,'ringLevels',0));
g = stagnant_cap_geometry(o);
s = stagnant_cap_system(g,o,struct('kind','force','Fg',[0 0 1]));
old = drop_recompleted_matrix(g.innerSources,g.outerSources, ...
    g.points,g.normals,o.muIn,o.completionScale);
assert(norm(s.matrix-old,'fro')<1e-10);
assert(norm(s.KN.'*s.lambda0-[0;0;1])<1e-14);
o = stagnant_cap_options(struct('sourceRequest',50,'phi',pi/2, ...
    'ringLevels',2));
g = stagnant_cap_geometry(o);
s = stagnant_cap_system(g,o,struct('kind','force','Fg',[0 0 1]));
v = [1;-2;3];
response = s.matrix*[s.KN*v;zeros(3*g.N,1)];
M = g.collocationCount;
Mc = nnz(g.cap);
assert(norm(response(1:3*M),inf)<1e-12);
expectedNormal = o.completionScale*g.N*g.normals*v;
assert(norm(response(3*M+1:4*M)-expectedNormal,inf)<1e-11);
capV = repmat(v.',Mc,1);
expectedTangential = [sum(capV.*s.tangentOne(g.cap,:),2), ...
    sum(capV.*s.tangentTwo(g.cap,:),2)];
expectedTangential = reshape(expectedTangential.',[],1);
assert(norm(response(4*M+1:4*M+2*Mc)- ...
    o.completionScale*g.N*expectedTangential,inf)<1e-11);
assert(norm(response(4*M+2*Mc+1:end),inf)<1e-12);
% The projected and prescribed-speed matrices must give the same physical
% boundary residual for matching lambda and U, including nonzero lambda0.
N = g.N;
alpha = sin((1:3*N).');
alphaNodes = reshape(alpha,3,N);
alphaNodes = alphaNodes-mean(alphaNodes,2)+repmat([0;0;0.2],1,N);
alpha = alphaNodes(:);
gamma = cos((1:3*N).');
lambdaNodes = alphaNodes-mean(alphaNodes,2)+repmat([0;0;1/N],1,N);
lambda = lambdaNodes(:);
Uz = -o.completionScale*sum(alphaNodes(3,:));
speed = stagnant_cap_system(g,o,struct('kind','speed','Uz',Uz));
projectedResidual = s.matrix*[alpha;gamma]-s.rhs;
physicalResidual = speed.matrix*[lambda;gamma]-speed.rhs;
assert(norm(projectedResidual-physicalResidual,inf)<1e-11);
% Added S and R blocks contribute force and torque with the same ordering
% used by the image-source kernels.
o.refinementFamilies = ["S","R"];
o.n_level = 1;
o.n_az_source = 8;
o.n_az_colloc = 12;
refined = stagnant_cap_refined_geometry(o);
mixed = stagnant_cap_system(refined,o, ...
    struct('kind','force','Fg',[0 0 1]));
assert(isequal(size(mixed.matrix),[6*refined.collocationCount, ...
    refined.exteriorCoefficientCount+refined.interiorCoefficientCount]));
z = zeros(mixed.outerColumns,1);
sourceStart = 3*refined.N;
rotletStart = sourceStart+3*o.n_az_source;
z(sourceStart+1) = 2;
z(rotletStart+2) = 3;
assert(norm(mixed.forceMap*z-[2;0;0])<1e-13);
expectedTorque = cross(refined.innerAddedSources(1,:)-o.centre, ...
    [2 0 0]).'+[0;3;0];
assert(norm(mixed.momentMap(4:6,:)*z-expectedTorque)<1e-12);
assert(norm(mixed.forceMap*mixed.lambda0-[0;0;1])<1e-13);
alphaMixed = sin((1:mixed.outerColumns).');
alphaMixed = alphaMixed- ...
    mixed.forceCompletion*(mixed.forceMap*alphaMixed)+ ...
    mixed.forceCompletion*[0;0;2];
gammaMixed = cos((1:mixed.innerColumns).');
lambdaMixed = alphaMixed- ...
    mixed.forceCompletion*(mixed.forceMap*alphaMixed)+mixed.lambda0;
UzMixed = -o.completionScale*mixed.forceMap*alphaMixed;
assert(norm(UzMixed(1:2))<1e-12);
direct = stagnant_cap_system(refined,o, ...
    struct('kind','speed','Uz',UzMixed(3)));
assert(norm(mixed.matrix*[alphaMixed;gammaMixed]-mixed.rhs- ...
    (direct.matrix*[lambdaMixed;gammaMixed]-direct.rhs),inf)<1e-10);
o.refinementFamilies = "S";
o.derivative_direction = "cone_normal";
derivativeGeometry = stagnant_cap_refined_geometry(o);
derivativeSystem = stagnant_cap_system(derivativeGeometry,o, ...
    struct('kind','force','Fg',[0 0 1]));
derivativeStart = 3*(derivativeGeometry.N+o.n_az_source);
probe = zeros(derivativeSystem.outerColumns,1);
probe(derivativeStart+3) = 2;
q = derivativeGeometry.innerDerivativeDirections(1,:);
d = derivativeGeometry.innerDerivativeDistances(1);
assert(norm(derivativeSystem.forceMap*probe)<1e-13);
assert(norm(derivativeSystem.momentMap(4:6,:)*probe- ...
    d*cross(q,[0 0 2]).')<1e-13);
fprintf('stagnant_cap_system self-test passed.\n');
end
