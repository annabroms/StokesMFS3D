function A = drop_image_interaction_block(source,target,t1,t2,mu)
%DROP_IMAGE_INTERACTION_BLOCK Incident normal velocity and traction jump.
if nargin == 0
    self_test(); A = []; return
end
x = target.collocation; n = target.normals; M = size(x,1);
be = source.basis;
ne = 3*sum(arrayfun(@(b) size(b.points,1),be));
ni = 3*sum(arrayfun(@(b) size(b.points,1),source.interiorBasis));
H = zeros(3,ne); offset = 0;
for k = 1:numel(be)
    count = size(be(k).points,1);
    if be(k).family=="S"
        H(:,offset+(1:3*count)) = be(k).scale*repmat(eye(3),1,count);
    end
    offset = offset+3*count;
end
C = H'/(H*H');
G = stokes_mfs_velocity_matrix(be,x);
T = stokes_mfs_traction_matrix(be,x,n);
Gp = G-(G*C)*H;
Tp = T-(T*C)*H;
A = zeros(6*M,ne+ni);
for m = 1:M
    rows = 3*m-2:3*m;
    A(3*M+m,1:ne) = n(m,:)*Gp(rows,:);
    A(4*M+2*m-1:4*M+2*m,1:ne) = ...
        (1-mu)*[t1(m,:);t2(m,:)]*Tp(rows,:);
end
end

function self_test()
% Match the existing pure-Stokeslet cross block.
v = [eye(3);-eye(3)]; x = [1 1 1;1 -1 -1;-1 1 -1;-1 -1 1]/sqrt(3);
b.family = "S"; b.points = [3 0 0]+0.6*v;
b.normalDirections = zeros(0,3); b.scale = 1;
s.basis = b; b.points = [3 0 0]+2*v; s.interiorBasis = b;
t.collocation = x; t.normals = x;
r = repmat([0 0 1],4,1);
t1 = cross(x,r,2); t1 = t1./vecnorm(t1,2,2); t2 = cross(x,t1,2);
A = drop_image_interaction_block(s,t,t1,t2,2);
B = drop_interaction_block(s.basis.points,x,x,t1,t2,2);
assert(norm(A-B,'fro')<1e-12);
fprintf('drop_image_interaction_block self-test passed.\n');
end
