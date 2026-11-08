function operators = bubble_boundary_operators(geometry,opt)
%BUBBLE_BOUNDARY_OPERATORS Shared unprojected exterior kernels and local frames.
% Cartesian rows/columns are node-interleaved; tangent rows are [t1_i;t2_i].
if nargin==0
    self_test();
    operators = [];
    return
end
M = geometry.collocationCount;
C = geometry.coefficientCount;
cap = find(geometry.cap);
mobile = find(geometry.mobile);
[basis,count] = bubble_basis(geometry,opt);
assert(count==C && numel(cap)+numel(mobile)==M);
G = stokes_mfs_velocity_matrix(basis,geometry.points)/opt.muOut;
Tm = stokes_mfs_traction_matrix(basis,geometry.points(mobile,:), ...
    geometry.normals(mobile,:));
assert(isequal(size(G),[3*M,C]) && isequal(size(Tm),[3*numel(mobile),C]));

% Preserve the previous solver's reference-axis choice and tangent ordering.
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
operators.basis = basis;
operators.G = G;
operators.Gc = G(capVectorRows,:);
operators.mobileVectorRows = mobileVectorRows;
operators.Tm = Tm;
operators.normalMap = normalMap;
operators.Nm = normalMap(mobile,mobileVectorRows);
operators.Dc = tangentMap(capTangentRows,capVectorRows);
operators.Dm = tangentMap(mobileTangentRows,mobileVectorRows);
operators.KM = kron(ones(M,1),eye(3));
operators.KMc = kron(ones(numel(cap),1),eye(3));
operators.KMm = kron(ones(numel(mobile),1),eye(3));
operators.cap = cap;
operators.mobile = mobile;
operators.tangentOne = t1;
operators.tangentTwo = t2;
end

function self_test()
% Verify kernel layouts, local orthonormal frames and both source-site modes.
for levels = [0 1]
    o = bubble_options(struct('sourceRequest',12,'phi',pi/3, ...
        'n_level',levels,'n_az_source',6,'n_az_colloc',8));
    g = bubble_geometry(o);
    p = bubble_boundary_operators(g,o);
    assert(isequal(size(p.G),[3*g.collocationCount,g.coefficientCount]));
    assert(norm(p.Nm*p.KMm-g.normals(g.mobile,:),'fro')<1e-13);
    for k = 1:g.collocationCount
        Q = [g.normals(k,:);p.tangentOne(k,:);p.tangentTwo(k,:)];
        assert(norm(Q*Q'-eye(3),'fro')<1e-12);
    end
    assert(norm(p.G-stokes_mfs_velocity_matrix(p.basis,g.points)/o.muOut,'fro')==0);
end
fprintf('bubble_boundary_operators self-test passed.\n');
end
