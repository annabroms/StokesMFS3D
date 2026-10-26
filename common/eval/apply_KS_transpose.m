function [Sbody,s5body,Ssum,s5sum] = apply_KS_transpose(rvec_in,q,lambda)
%APPLY_KS_TRANSPOSE Extract force-on-fluid stresslets from Stokeslet strengths.
%   [Sbody,s5body,Ssum,s5sum] = apply_KS_transpose(rvec_in,q,lambda)
% rvec_in is (N*P)-by-3, q is P-by-3, and lambda is interleaved
% [fx1;fy1;fz1;fx2;...] by body. Lambda contains point-force amplitudes,
% not surface force densities, so no quadrature weights are applied.
% Sbody is 3-by-3-by-P; s5body is 5-by-P in the orthonormal Frobenius basis
% E1=diag(1,-1,0)/sqrt(2), E2=diag(1,1,-2)/sqrt(6),
% E3=(ex*ey'+ey*ex')/sqrt(2), E4=(ex*ez'+ez*ex')/sqrt(2),
% E5=(ey*ez'+ez*ey')/sqrt(2). Ssum and s5sum sum over bodies.

if nargin == 0
    self_test();
    Sbody = []; s5body = []; Ssum = []; s5sum = [];
    return
end
validateattributes(rvec_in,{'numeric'},{'real','finite','2d','ncols',3});
validateattributes(q,{'numeric'},{'real','finite','2d','ncols',3});
P = size(q,1);
assert(P >= 1 && mod(size(rvec_in,1),P)==0, ...
    'apply_KS_transpose:BadBodyLayout');
N = size(rvec_in,1)/P;
validateattributes(lambda,{'numeric'}, ...
    {'real','finite','column','numel',3*N*P});

E = zeros(3,3,5);
E(:,:,1) = diag([1,-1,0])/sqrt(2);
E(:,:,2) = diag([1,1,-2])/sqrt(6);
E(1,2,3) = 1/sqrt(2); E(2,1,3) = E(1,2,3);
E(1,3,4) = 1/sqrt(2); E(3,1,4) = E(1,3,4);
E(2,3,5) = 1/sqrt(2); E(3,2,5) = E(2,3,5);

Sbody = zeros(3,3,P);
s5body = zeros(5,P);
for p = 1:P
    nodes = (p-1)*N+(1:N);
    rho = rvec_in(nodes,:)-q(p,:);
    force = reshape(lambda(3*(nodes(1)-1)+1:3*nodes(end)),3,[]).';
    moment = rho.'*force;
    Sbody(:,:,p) = (moment+moment.')/2-trace(moment)*eye(3)/3;
    for alpha = 1:5
        s5body(alpha,p) = sum(Sbody(:,:,p).*E(:,:,alpha),'all');
    end
end
Ssum = sum(Sbody,3);
s5sum = sum(s5body,2);
end

function self_test()
%SELF_TEST checks force moments and their far-field Stokeslet dipole.
q = [0.3,-0.2,0.1;-0.4,0.5,-0.3];
rho = [0.4,0.1,-0.2;-0.3,0.2,0.1];
r = [q(1,:)+rho;q(2,:)+rho];
lambda = [1;2;3;-1;0;2;2;-1;1;0;1;-2];
[S,v,St,vt] = apply_KS_transpose(r,q,lambda);
assert(isequal(size(S),[3,3,2]) && isequal(size(v),[5,2]));
assert(norm(St-sum(S,3),'fro')<1e-14 && norm(vt-sum(v,2))<1e-14);
for p = 1:2
    assert(norm(S(:,:,p)-S(:,:,p).','fro')<1e-14);
    assert(abs(trace(S(:,:,p)))<1e-14);
end
[Slinear,~,~,~] = apply_KS_transpose(r,q,2*lambda);
assert(norm(Slinear(:)-2*S(:))<1e-13);
shift = [0.7,-1.1,0.2];
[Sshift,~,~,~] = apply_KS_transpose(r+shift,q+shift,lambda);
assert(norm(Sshift(:)-S(:))<1e-13);
angle = 0.41;
Q = [cos(angle),-sin(angle),0;sin(angle),cos(angle),0;0,0,1];
rotatedR = (Q*r.').';
rotatedQ = (Q*q.').';
rotatedLambda = reshape(Q*reshape(lambda,3,[]),[],1);
[Srot,~,~,~] = apply_KS_transpose(rotatedR,rotatedQ,rotatedLambda);
for p = 1:2
    assert(norm(Srot(:,:,p)-Q*S(:,:,p)*Q.','fro')<1e-13);
end

% A force-neutral pair has a known leading exterior-flow stresslet.
d = 0.4;
src = [d/2,0,0;-d/2,0,0];
f = [1;0;0;-1;0;0];
[Sdip,s5dip,~,~] = apply_KS_transpose(src,[0,0,0],f);
x = [100,17,-13;55,-130,21;-42,29,137; ...
    -111,-53,47;67,99,-81;83,-74,-112; ...
    -96,113,31;121,51,79];
u = reshape(stokes_SLP_mat(src,x)*f,3,[]).';
uFar = zeros(size(x));
E = zeros(3,3,5);
E(:,:,1) = diag([1,-1,0])/sqrt(2);
E(:,:,2) = diag([1,1,-2])/sqrt(6);
E(1,2,3) = 1/sqrt(2); E(2,1,3) = E(1,2,3);
E(1,3,4) = 1/sqrt(2); E(3,1,4) = E(1,3,4);
E(2,3,5) = 1/sqrt(2); E(3,2,5) = E(2,3,5);
fitMatrix = zeros(3*size(x,1),5);
for j = 1:size(x,1)
    xx = x(j,:).';
    uFar(j,:) = (3/(8*pi))* ...
        (xx*(xx.'*Sdip(:,:,1)*xx)/norm(xx)^5).';
    for alpha = 1:5
        fitMatrix(3*(j-1)+(1:3),alpha) = ...
            (3/(8*pi))*xx*(xx.'*E(:,:,alpha)*xx)/norm(xx)^5;
    end
end
assert(norm(u-uFar,'fro')/norm(uFar,'fro')<1e-4);
fittedStresslet = fitMatrix\reshape(u.',[],1);
assert(norm(fittedStresslet-s5dip)/norm(s5dip)<1e-4);
fprintf('apply_KS_transpose self-test passed.\n');
end
