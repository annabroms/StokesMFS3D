function system = bubble_system(geometry,opt,Fg)
%BUBBLE_SYSTEM Assemble the projected exterior mixed-boundary system.
% Preserve interleaved [x1;y1;z1;x2;...] source and target ordering.
if nargin == 0
    bubble_self_test();
    system = [];
    return
end
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
[basis,coefficientCount] = bubble_basis(geometry,opt);
assert(coefficientCount==N);
G = stokes_mfs_velocity_matrix(basis,geometry.points)/opt.muOut;
Tm = stokes_mfs_traction_matrix(basis,geometry.points(mobile,:), ...
    geometry.normals(mobile,:));
assert(isequal(size(G),[3*M,N]) && isequal(size(Tm),[3*Mm,N]));

reference = repmat([1 0 0],M,1);
reference(abs(geometry.normals(:,1))>0.9,:) = ...
    repmat([0 1 0],nnz(abs(geometry.normals(:,1))>0.9),1);
t1 = cross(geometry.normals,reference,2);
t1 = t1./vecnorm(t1,2,2);
t2 = cross(geometry.normals,t1,2);
normalMap = sparse(repelem((1:M).',3), ...
    (1:3*M).',reshape(geometry.normals.',[],1),M,3*M);
tangentRows = zeros(6*M,1);
tangentColumns = zeros(6*M,1);
tangentValues = zeros(6*M,1);
for p = 1:M
    rows = 6*p-5:6*p;
    columns = (3*p-2:3*p).';
    tangentRows(rows) = [repmat(2*p-1,3,1);repmat(2*p,3,1)];
    tangentColumns(rows) = [columns;columns];
    tangentValues(rows) = [t1(p,:).';t2(p,:).'];
end
tangentMap = sparse(tangentRows,tangentColumns,tangentValues,2*M,3*M);
capVectorRows = reshape([3*cap-2,3*cap-1,3*cap].',[],1);
mobileVectorRows = reshape([3*mobile-2,3*mobile-1,3*mobile].',[],1);
capTangentRows = reshape([2*cap-1,2*cap].',[],1);
mobileTangentRows = reshape([2*mobile-1,2*mobile].',[],1);
Dc = tangentMap(capTangentRows,capVectorRows);
Dm = tangentMap(mobileTangentRows,mobileVectorRows);
Gc = G(capVectorRows,:);
KM = kron(ones(M,1),eye(3));
KMc = kron(ones(Mc,1),eye(3));

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
system.matrix = [normalMap*Gp+beta*(normalMap*KM)*K; ...
    Dc*Gcp+beta*(Dc*KMc)*K; Dm*Tmp];
system.rhs = -[normalMap*G*lambda0; ...
    Dc*Gc*lambda0; Dm*Tm*lambda0];
system.forceMap = K;
system.forceCompletion = C;
system.lambda0 = lambda0;
system.momentMap = moments;
system.basis = basis;
system.normalRows = M;
system.capTangentRows = 2*Mc;
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
