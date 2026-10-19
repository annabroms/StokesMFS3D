function [A,H,t1,t2] = drop_image_self_block(body,mu,c)
%DROP_IMAGE_SELF_BLOCK Assemble one drop's mixed-source boundary equations.
% Exterior strengths are projected to zero net Stokeslet force; interior
% source strengths remain unprojected. Rows retain [continuity;normal;tangent].
if nargin == 0
    self_test(); A = []; H = []; t1 = []; t2 = []; return
end
x = body.collocation;
n = body.normals;
M = size(x,1);
R = repmat([1 0 0],M,1);
R(abs(n(:,1))>0.9,:) = repmat([0 1 0],sum(abs(n(:,1))>0.9),1);
t1 = cross(n,R,2); t1 = t1./vecnorm(t1,2,2);
t2 = cross(n,t1,2);
be = body.basis;
bi = body.interiorBasis;
ne = 3*sum(arrayfun(@(b) size(b.points,1),be));
ni = 3*sum(arrayfun(@(b) size(b.points,1),bi));
H = zeros(3,ne);
offset = 0;
for k = 1:numel(be)
    count = size(be(k).points,1);
    if be(k).family=="S"
        H(:,offset+(1:3*count)) = ...
            be(k).scale*repmat(eye(3),1,count);
    end
    offset = offset+3*count;
end
C = H'/(H*H');
G = stokes_mfs_velocity_matrix(be,x);
E = stokes_mfs_velocity_matrix(bi,x)/mu;
T = stokes_mfs_traction_matrix(be,x,n);
Q = stokes_mfs_traction_matrix(bi,x,n);
Gp = G-(G*C)*H;
Tp = T-(T*C)*H;
A = zeros(6*M,ne+ni);
A(1:3*M,:) = [Gp,-E];
for m = 1:M
    rows = 3*m-2:3*m;
    A(3*M+m,1:ne) = n(m,:)*Gp(rows,:)+c*n(m,:)*H;
    A(4*M+2*m-1:4*M+2*m,1:ne) = ...
        [t1(m,:);t2(m,:)]*Tp(rows,:);
    A(4*M+2*m-1:4*M+2*m,ne+(1:ni)) = ...
        -[t1(m,:);t2(m,:)]*Q(rows,:);
end
end

function self_test()
% Recover the basic Stokeslet block and check all-image force projection.
v = [eye(3);-eye(3)]; x = [1 1 1;1 -1 -1;-1 1 -1;-1 -1 1]/sqrt(3);
b.family = "S"; b.points = 0.6*v; b.normalDirections = zeros(0,3); b.scale = 1;
body.basis = b; b.points = 2*v; body.interiorBasis = b;
body.collocation = x; body.normals = x;
[A,H] = drop_image_self_block(body,2,1/(8*pi));
[B,~,~] = drop_recompleted_matrix(0.6*v,2*v,x,x,2,1/(8*pi));
assert(norm(A-B,'fro')<1e-12 && isequal(size(H),[3 18]));
b.points = [0.8 0 0]; body.basis(2) = b;
[A,H] = drop_image_self_block(body,2,1/(8*pi));
assert(isequal(size(A),[24 39]));
C = H'/(H*H');
z = sin((1:size(H,2)).');
assert(norm(H*(z-C*(H*z)))<1e-12);
body.basis(2).scale = 2;
for k = 3:5
    body.basis(k) = body.basis(2);
end
body.basis(3).family = "R";
body.basis(4).family = "T";
body.basis(4).normalDirections = [1 0 0];
body.basis(5).family = "D";
[~,H] = drop_image_self_block(body,2,1/(8*pi));
assert(norm(H(:,19:21)-2*eye(3),'fro')<1e-14);
assert(~any(H(:,22:end),'all'));
body.basis(4).family = "T_n";
body.basis(6) = body.basis(4);
body.basis(6).family = "T_t";
body.basis(6).normalDirections = [0 1 0];
[~,H] = drop_image_self_block(body,2,1/(8*pi));
if any(H(:,22:end),'all')
    error('drop_image_self_block:SelfTest','Only Stokeslets may enter the force map.');
end
fprintf('drop_image_self_block self-test passed.\n');
end
