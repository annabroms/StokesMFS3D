function system = bubble_system(geometry,opt,Fg,capMode)
%BUBBLE_SYSTEM Assemble the projected exterior mixed-boundary system.
% Preserve interleaved [x1;y1;z1;x2;...] source and target ordering.
if nargin == 0
    bubble_self_test();
    system = [];
    return
end
if nargin<4, capMode = "split"; end
capMode = string(capMode);
assert(isscalar(capMode) && ismember(capMode,["split","cartesian"]));
Fg = Fg(:);
assert(numel(Fg)==3 && all(isfinite(Fg)) && isreal(Fg) && ...
    norm(Fg)>0 && norm(Fg(1:2))<1e-13*norm(Fg), ...
    'Prescribed force must be nonzero and axial.');
M = geometry.collocationCount;
N = geometry.coefficientCount;
cap = find(geometry.cap);
mobile = find(geometry.mobile);
Mc = numel(cap);
Mm = numel(mobile);
operators = bubble_boundary_operators(geometry,opt);
basis = operators.basis;
G = operators.G;
Tm = operators.Tm;
normalMap = operators.normalMap;
Dc = operators.Dc;
Dm = operators.Dm;
Gc = operators.Gc;
KM = operators.KM;
KMc = operators.KMc;

body.basis = basis;
[moments,~] = image_mobility_maps(body,opt.centre);
K = moments(1:3,:);
C = K'/(K*K');
assert(isequal(size(K),[3,N]) && norm(K*C-eye(3),'fro')<1e-11);
lambda0 = C*Fg;
Gp = G-(G*C)*K;
Gcp = Gc-(Gc*C)*K;
Tmp = Tm-(Tm*C)*K;
beta = opt.completionScale;
if capMode=="split"
    system.matrix = [normalMap*Gp+beta*(normalMap*KM)*K; ...
        Dc*Gcp+beta*(Dc*KMc)*K; Dm*Tmp];
    system.rhs = -[normalMap*G*lambda0; ...
        Dc*Gc*lambda0; Dm*Tm*lambda0];
    system.normalRows = M;
    system.capTangentRows = 2*Mc;
    system.capVelocityRows = 0;
    system.blockCounts = [M,2*Mc,2*Mm];
else
    % All physical unknown source forces remain projected; only cap rows change.
    Bm = operators.Nm*operators.KMm;
    system.matrix = [Gcp+beta*KMc*K; ...
        operators.Nm*Gp(operators.mobileVectorRows,:)+beta*Bm*K; Dm*Tmp];
    system.rhs = -[Gc*lambda0; ...
        operators.Nm*G(operators.mobileVectorRows,:)*lambda0; Dm*Tm*lambda0];
    system.normalRows = Mm;
    system.capTangentRows = 0;
    system.capVelocityRows = 3*Mc;
    system.blockCounts = [3*Mc,Mm,2*Mm];
end
system.capMode = capMode;
system.forceMap = K;
system.forceCompletion = C;
system.lambda0 = lambda0;
system.momentMap = moments;
system.basis = basis;
system.mobileTractionRows = 2*Mm;
system.sourceCount = geometry.sourceCount;
system.coefficientCount = N;
system.collocationCount = M;
assert(isequal(size(system.matrix),[3*M,N]) && ...
    isequal(size(system.rhs),[3*M,1]));
end

function bubble_self_test()
% Compare the projected rows with explicit force and translation substitution.
for families = ["S","SRT_t"]
    o = bubble_options(struct('sourceRequest',50,'phi',pi/3, ...
        'n_level',1,'n_az_source',6,'n_az_colloc',8, ...
        'refinementFamilies',families));
    g = bubble_geometry(o);
    Fg = [0;0;-1];
    s = bubble_system(g,o,Fg);
    assert(isequal(size(s.matrix),[3*g.collocationCount, ...
        g.coefficientCount]));
    assert(norm(s.forceMap*s.lambda0-Fg)<1e-13);
    assert(any(s.forceMap(3,3*g.baseSourceCount+1: ...
        3*(g.baseSourceCount+g.coneSourceCount))));
    alpha = sin((1:g.coefficientCount).');
    lambda = alpha-s.forceCompletion*(s.forceMap*alpha)+s.lambda0;
    U = -o.completionScale*s.forceMap*alpha;
    basis = s.basis;
    G = stokes_mfs_velocity_matrix(basis,g.points)/o.muOut;
    T = stokes_mfs_traction_matrix(basis,g.points(g.mobile,:), ...
        g.normals(g.mobile,:));
    velocity = reshape(G*lambda,3,[]).'-U.';
    normal = sum(velocity.*g.normals,2);
    capVelocity = velocity(g.cap,:);
    capNormal = g.normals(g.cap,:);
    capTangent = capVelocity-sum(capVelocity.*capNormal,2).*capNormal;
    traction = reshape(T*lambda,3,[]).';
    mobileNormal = g.normals(g.mobile,:);
    mobileTangent = traction-sum(traction.*mobileNormal,2).*mobileNormal;
    expectedNorm = norm([normal;capTangent(:);mobileTangent(:)]);
    actualNorm = norm(s.matrix*alpha-s.rhs);
    assert(abs(expectedNorm-actualNorm)<1e-10*max(1,actualNorm));
    cartesian = bubble_system(g,o,Fg,"cartesian");
    cartResidual = cartesian.matrix*alpha-cartesian.rhs;
    capCartesian = reshape(velocity(g.cap,:).',[],1);
    assert(norm(cartResidual(1:3*g.capCount)-capCartesian)<1e-10);
    assert(norm(cartResidual(3*g.capCount+(1:g.mobileCount))-normal(g.mobile))<1e-10);
    assert(abs(norm(cartResidual)-actualNorm)<1e-10*max(1,actualNorm));
    assert(isequal(cartesian.blockCounts,[3*g.capCount,g.mobileCount,2*g.mobileCount]));
    assert(cartesian.capTangentRows==0 && cartesian.normalRows==g.mobileCount);
end
derivativeOptions = bubble_options(struct('sourceRequest',50, ...
    'phi',pi/3,'n_level',1,'n_az_source',6,'n_az_colloc',8, ...
    'refinementFamilies',["S","D","Q_s"], ...
    'derivative_direction',"radial"));
derivativeGrid = bubble_geometry(derivativeOptions);
derivativeSystem = bubble_system( ...
    derivativeGrid,derivativeOptions,[0;0;-1]);
assert(isequal(size(derivativeSystem.matrix), ...
    [3*derivativeGrid.collocationCount, ...
    derivativeGrid.coefficientCount]));
nonStokesletStart = 3*(derivativeGrid.baseSourceCount+ ...
    derivativeGrid.coneSourceCount)+1;
assert(all(derivativeSystem.forceMap(:,nonStokesletStart:end)==0,'all'));
fprintf('bubble_system self-test passed.\n');
end
