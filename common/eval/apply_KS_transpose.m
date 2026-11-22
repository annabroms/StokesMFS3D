function [Sbody,s5body,Ssum,s5sum,parts] = apply_KS_transpose(rvec_in,q,lambda)
%APPLY_KS_TRANSPOSE Compute each body's force-on-fluid stresslet K_S^T lambda.
%
% Ordinary Stokeslet sources:
%   [Sbody,s5body,Ssum,s5sum,parts] = apply_KS_transpose(rvec_in,q,lambda)
%
%   Inputs (P bodies, N sources per body):
%     rvec_in  (N*P)-by-3 source positions, grouped by body.
%     q        P-by-3 body centres; moments are taken about these centres.
%     lambda   (3*N*P)-by-1 physical point-force amplitudes, grouped by body
%              and interleaved [fx1;fy1;fz1;fx2;fy2;fz2;...].
%
% Image-enriched sources:
%   [Sbody,s5body,Ssum,s5sum,parts] = ...
%       apply_KS_transpose(details.grid,details.lambda)
%
%   Inputs (from an image mobility or resistance solve):
%     grid     grid.centres is P-by-3. grid.body(p).basis contains the
%              source blocks, with fields family, points, normalDirections
%              (for T blocks), and scale. Source counts may differ by body.
%     lambda   P-element cell array of physical source-strength columns,
%              ordered by basis block, then interleaved xyz within a block.
%              Use the solver's projected and recompleted strengths.
%              In this two-input call, lambda occupies the argument q.
%
% Outputs (identical for both input forms):
%   Sbody      3-by-3-by-P symmetric, trace-free stresslet tensors.
%   s5body     5-by-P coordinates: s5body(alpha,p) = <Sbody(:,:,p),E_alpha>_F.
%   Ssum       3-by-3 tensor, summed over all bodies.
%   s5sum      5-by-1 coordinates, summed over all bodies.
%   parts      Decomposition into two 3-by-3-by-P arrays:
%                .stokeslet          moments of S sources about each centre;
%                .intrinsicStresslet contributions from T sources.
%              Sbody = parts.stokeslet + parts.intrinsicStresslet.
%
% Conventions:
%   STF(A) = (A+A.')/2 - trace(A)*I/3.
%   An S source at y with strength c contributes STF((y-q)*c.').
%   A T source with direction n contributes -2*STF(n*c.').
%   R and D sources contribute zero to this symmetric force-dipole moment.
%   Each block's scale is applied once by image_source_coefficients.
%   Strengths are amplitudes, so no surface quadrature weights are applied.
%   These are force-on-fluid moments, before any suspension-stress sign,
%   volume normalization, or thermal-energy factor is applied.
%
% The orthonormal Frobenius basis uses Cartesian unit vectors ex, ey, ez:
%   E1 = diag(1,-1,0)/sqrt(2), E2 = diag(1,1,-2)/sqrt(6),
%   E3 = (ex*ey.'+ey*ex.')/sqrt(2),
%   E4 = (ex*ez.'+ez*ex.')/sqrt(2),
%   E5 = (ey*ez.'+ez*ey.')/sqrt(2).
%
% Call apply_KS_transpose() to run the small deterministic self-test.

%% Self-test entry point
if nargin == 0
    self_test();
    Sbody = [];
    s5body = [];
    Ssum = [];
    s5sum = [];
    parts = [];
    return
end

%% Read the input form and check the body/source layout
imageForm = isstruct(rvec_in);
if imageForm
    assert(nargin == 2 && isfield(rvec_in,'body') && ...
        isfield(rvec_in,'centres') && iscell(q), ...
        'Image input must be (details.grid,details.lambda).');
    grid = rvec_in;
    strengths = q;
    q = grid.centres;
else
    assert(nargin == 3,'Expected source nodes, centres, and strengths.');
    validateattributes(rvec_in,{'numeric'},{'real','finite','2d','ncols',3});
end
validateattributes(q,{'numeric'},{'real','finite','2d','ncols',3});
P = size(q,1);
if imageForm
    assert(P >= 1 && numel(grid.body)==P && numel(strengths)==P, ...
        'apply_KS_transpose:BadBodyLayout');
else
    assert(P >= 1 && mod(size(rvec_in,1),P)==0, ...
        'apply_KS_transpose:BadBodyLayout');
    N = size(rvec_in,1)/P;
    validateattributes(lambda,{'numeric'}, ...
        {'real','finite','column','numel',3*N*P});
end

%% Accumulate physical source moments about each body centre
% Keep the two source representations separate; both fill the same parts.
parts = struct('stokeslet',zeros(3,3,P), ...
    'intrinsicStresslet',zeros(3,3,P));
if imageForm
    for p = 1:P
        basis = grid.body(p).basis;
        coefficientCount = 3*sum(arrayfun(@(b) size(b.points,1),basis));
        validateattributes(strengths{p},{'numeric'}, ...
            {'real','finite','column','numel',coefficientCount});
        coefficients = image_source_coefficients(basis,strengths{p});

        for b = 1:numel(basis)
            block = basis(b);
            blockStrengths = coefficients{b};
            if block.family == "S"
                sourceOffsets = block.points-q(p,:);
                moment = sourceOffsets.'*blockStrengths;
                parts.stokeslet(:,:,p) = parts.stokeslet(:,:,p)+ ...
                    stf(moment);

            elseif block.family == "T"
                assert(isequal(size(block.normalDirections),size(blockStrengths)));
                moment = block.normalDirections.'*blockStrengths;
                parts.intrinsicStresslet(:,:,p) = ...
                    parts.intrinsicStresslet(:,:,p)-2*stf(moment);

            else
                % Rotlets and source dipoles have no STF force-dipole moment.
                assert(block.family == "R" || block.family == "D", ...
                    'apply_KS_transpose:UnknownFamily');
            end
        end
    end
else
    for p = 1:P
        sourceNodes = (p-1)*N+(1:N);
        sourceRows = 3*(sourceNodes(1)-1)+1:3*sourceNodes(end);
        sourceOffsets = rvec_in(sourceNodes,:)-q(p,:);
        forces = reshape(lambda(sourceRows),3,[]).';
        parts.stokeslet(:,:,p) = stf(sourceOffsets.'*forces);
    end
end
Sbody = parts.stokeslet+parts.intrinsicStresslet;

%% Express the tensors in the five-component orthonormal basis
E = zeros(3,3,5);
E(:,:,1) = diag([1,-1,0])/sqrt(2);
E(:,:,2) = diag([1,1,-2])/sqrt(6);
E(1,2,3) = 1/sqrt(2);
E(2,1,3) = E(1,2,3);
E(1,3,4) = 1/sqrt(2);
E(3,1,4) = E(1,3,4);
E(2,3,5) = 1/sqrt(2);
E(3,2,5) = E(2,3,5);

s5body = zeros(5,P);
for p = 1:P
    for alpha = 1:5
        s5body(alpha,p) = sum(Sbody(:,:,p).*E(:,:,alpha),'all');
    end
end

%% Sum over bodies
Ssum = sum(Sbody,3);
s5sum = sum(s5body,2);
end

function S = stf(moment)
%STF removes the antisymmetric and isotropic parts of a force moment.
S = (moment+moment.')/2-trace(moment)*eye(3)/3;
end

function self_test()
%SELF_TEST checks force moments and their far-field Stokeslet dipole.

%% Ordinary sources: layout, symmetry, linearity and rigid transformations
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

%% Ordinary sources: compare the stresslet with its far-field velocity
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

%% Image sources: intrinsic T moment and scaled S/T decomposition
% The image T kernel also has an isotropic radial source-flow component;
% fit it separately from its five STF force-dipole coordinates.
block = struct('family',"T",'points',[0,0,0], ...
    'normalDirections',[1,0,0],'scale',1);
grid = struct('body',struct('basis',block), ...
    'centres',[0,0,0]);
strength = {[0.3;0.7;-0.2]};
[St,s5t,~,~,part] = apply_KS_transpose(grid,strength);
assert(norm(St(:,:,1)+2*stf( ...
    block.normalDirections.'*strength{1}.'))<1e-13);
assert(norm(part.intrinsicStresslet(:)-St(:))<1e-13);
assert(norm(part.stokeslet(:))==0);
stokesletBlock = block;
stokesletBlock.family = "S";
stokesletBlock.points = [0.2,0,0];
stokesletBlock.normalDirections = zeros(0,3);
stokesletBlock.scale = 2;
scaledT = block;
scaledT.scale = 0.5;
mixedGrid = struct('body',struct( ...
    'basis',[stokesletBlock,scaledT]),'centres',[0,0,0]);
mixedStrength = {[0;1;0;strength{1}]};
[mixedStress,~,~,~,mixedParts] = apply_KS_transpose( ...
    mixedGrid,mixedStrength);
expectedS = stf([0.2;0;0]*[0,2,0]);
expectedT = -stf(block.normalDirections.'*strength{1}.');
assert(norm(mixedParts.stokeslet(:,:,1)-expectedS,'fro')<1e-13);
assert(norm(mixedParts.intrinsicStresslet(:,:,1)-expectedT,'fro')<1e-13);
assert(norm(mixedStress(:,:,1)-expectedS-expectedT,'fro')<1e-13);
%% Image sources: recover the T moment from the velocity field
fitWithFlux = [fitMatrix,zeros(size(fitMatrix,1),1)];
for j = 1:size(x,1)
    xx = x(j,:).';
    fitWithFlux(3*(j-1)+(1:3),6) = xx/norm(xx)^3;
end
velocityT = stokes_mfs_velocity_matrix(block,x)*strength{1};
fitT = fitWithFlux\velocityT;
assert(norm(fitT(1:5)-s5t)/norm(s5t)<1e-11);
%% Image sources: rigid rotation and zero R/D contributions
G = [cos(0.37),-sin(0.37),0; ...
    sin(0.37),cos(0.37),0;0,0,1];
rotatedBlock = block;
rotatedBlock.points = (G*block.points.').';
rotatedBlock.normalDirections = ...
    (G*block.normalDirections.').';
rotatedGrid = struct('body',struct('basis',rotatedBlock), ...
    'centres',[0,0,0]);
rotatedStrength = {G*strength{1}};
[rotatedTensor,~,~,~] = apply_KS_transpose( ...
    rotatedGrid,rotatedStrength);
assert(norm(rotatedTensor(:,:,1)-G*St(:,:,1)*G.','fro')<1e-13);
for family = ["R","D"]
    block.family = family;
    [zeroStress,~,~,~] = apply_KS_transpose( ...
        struct('body',struct('basis',block),'centres',[0,0,0]), ...
        strength);
    assert(norm(zeroStress(:))==0);
end
fprintf('apply_KS_transpose self-test passed.\n');
end
