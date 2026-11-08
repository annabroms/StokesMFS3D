function [basis,C,w,k,moments,selected] = ring_mixed_basis(Y0,id0,Yc,idc,centre,radius,families,mode)
%RING_MIXED_BASIS Build axisymmetric mixed-source blocks and force completion.
% [basis,C,w,k,moments] = ring_mixed_basis(Y0,id0,Yc,idc,centre,radius,families)
% Y0,Yc are N0-by-3 and Nc-by-3 source coordinates, grouped by the
% corresponding contiguous 1-based ring indices id0,idc. centre is 1-by-3;
% radius is the physical radius; families are ring_family_list tokens.
% The GL block is always S; each selected family is placed at every clustered
% node. C maps ring-major [rho;psi;z] amplitudes per basis block to nodal
% Cartesian [x1;y1;z1;...] coefficients. w and k are p-by-1, w'*k=1;
% moments is 6-by-p, with force then torque about centre. No normalization
% is included in C. Q_s uses scale radius^2; other kernels use scale 1.
% T_t/T use the source ray, T_n the meridional normal, T_n2 the azimuthal
% normal. Unlike a fixed Cartesian frame these directions rotate with psi.
% Optional mode is "full" (default) or "swirl_free" ("swirl-free" also
% accepted for compatibility). In swirl-free mode,
% S, S_y, D, Q_s, T, T_t and T_n retain rho/z; R and T_n2 retain psi.
% S_y is d*(q dot grad_y)S, with d the interior source depth and q the
% meridional direction of increasing colatitude. Its coefficient has force
% units; it has no net force but may contribute torque.
% selected indexes these columns in the full block-major [rho;psi;z] map.
% The coefficient symmetry is exact; finite ring quadrature can leave small
% off-meridian azimuthal field components.
if nargin==0
    self_test();
    return
end
if nargin<8, mode = "full"; end
mode = string(mode);
if ~isscalar(mode) || ~ismember(mode,["full","swirl_free","swirl-free"])
    error('ring_mixed_basis:Mode','Use "full" or "swirl_free".');
end
if mode=="swirl-free", mode = "swirl_free"; end
if isempty(families)
    families = strings(1,0);
else
    families = ring_family_list(families);
end
[C0,~,~] = frame_map(Y0,id0,centre);
block = struct('family',"S",'points',Y0,'normalDirections',zeros(0,3), ...
    'scale',1,'derivativeDistances',zeros(0,1));
basis = block;
C = C0;
if ~isempty(families)
    [Cc,ray,meridian] = frame_map(Yc,idc,centre);
    psi = atan2(Yc(:,2)-centre(2),Yc(:,1)-centre(1));
    azimuth = [-sin(psi),cos(psi),zeros(numel(psi),1)];
    for family = families
        block.family = family;
        block.points = Yc;
        block.normalDirections = zeros(0,3);
        block.derivativeDistances = zeros(0,1);
        block.scale = 1;
        switch family
            case {"T","T_t"}
                block.normalDirections = ray;
            case "T_n"
                block.normalDirections = meridian;
            case "T_n2"
                block.normalDirections = azimuth;
            case "Q_s"
                block.normalDirections = ray;
                block.scale = radius^2;
            case "S_y"
                block.normalDirections = meridian;
                block.derivativeDistances = radius-vecnorm(Yc-centre,2,2);
                if any(block.derivativeDistances<=0)
                    error('ring_mixed_basis:SourceDepth', ...
                        'S_y sources must be strictly inside the sphere.');
                end
        end
        basis(end+1,1) = block; %#ok<AGROW>
        C = blkdiag(C,Cc);
    end
end
fullWidth = size(C,2);
selected = (1:fullWidth).';
if mode=="swirl_free"
    keep = false(size(selected));
    nr0 = max(id0);
    keep(1:3:3*nr0) = true;
    keep(3:3:3*nr0) = true;
    offset = 3*nr0;
    if ~isempty(families)
        nrc = max(idc);
        for family = families
            if ismember(family,["R","T_n2"])
                keep(offset+(2:3:3*nrc)) = true;
            else
                keep(offset+(1:3:3*nrc)) = true;
                keep(offset+(3:3:3*nrc)) = true;
            end
            offset = offset+3*nrc;
        end
    end
    selected = find(keep);
    C = C(:,selected);
end
body.basis = basis;
H = image_mobility_maps(body,centre);
moments = H*C;
w = moments(3,:).';
kFull = zeros(fullWidth,1);
kFull(3:3:3*max(id0)) = 1/size(Y0,1);
k = kFull(selected);
if abs(w.'*k-1)>1e-12
    error('ring_mixed_basis:Completion','Axial completion is not normalized.');
end
end

function [C,ray,meridian] = frame_map(Y,id,centre)
%FRAME_MAP Copy three rotating cylindrical amplitudes to each ring node.
id = id(:); N = size(Y,1); R = max(id);
if numel(id)~=N || ~isequal(unique(id),(1:R).')
    error('ring_mixed_basis:RingIndex','Ring indices must be contiguous and match the nodes.');
end
rel = Y-centre;
psi = atan2(rel(:,2),rel(:,1));
rho = hypot(rel(:,1),rel(:,2));
lengths = vecnorm(rel,2,2);
if any(rho==0 | lengths==0)
    error('ring_mixed_basis:AxisSource','Ring nodes must be off the axis.');
end
ray = rel./lengths;
meridian = [rel(:,3).*cos(psi)./lengths, ...
    rel(:,3).*sin(psi)./lengths,-rho./lengths];
j = (1:N).';
rows = [3*j-2;3*j-1;3*j-2;3*j-1;3*j];
cols = [3*id-2;3*id-2;3*id-1;3*id-1;3*id];
values = [cos(psi);sin(psi);-sin(psi);cos(psi);ones(N,1)];
C = sparse(rows,cols,values,3*N,3*R);
end

function self_test()
%SELF_TEST Verify copying, moments, completion and rotation of mixed kernels.
n = 8; psi = (0:n-1)'*2*pi/n;
Y = [0.4*cos(psi),0.4*sin(psi),0.2*ones(n,1)];
Yc = [0.6*cos(psi),0.6*sin(psi),0.4*ones(n,1)];
id = ones(n,1);
[basis,C,w,k,H] = ring_mixed_basis(Y,id,Yc,id,[0 0 0],1,["S","R","T_t","T_n","D","Q_s"]);
assert(isequal(size(C),[3*n*7,21]) && abs(w.'*k-1)<1e-14);
[original,~,~] = ring_cylindrical_maps(id,psi);
assert(norm(C(1:3*n,[1,3])-original,'fro')<1e-14);
[~,aliasMap] = ring_mixed_basis(Y,id,Yc,id,[0 0 0],1,"S","swirl-free");
[~,underscoredMap] = ring_mixed_basis(Y,id,Yc,id,[0 0 0],1,"S","swirl_free");
assert(isequal(aliasMap,underscoredMap));
v = zeros(21,1); v(1:3) = [2;3;4];
nodal = reshape(C*v,3,[]);
assert(norm(nodal(:,1)-[2;3;4])<1e-14 && norm(nodal(:,3)-[-3;2;4])<1e-14);
assert(norm(H(1:3,:)*k-[0;0;1])<1e-13);
v = zeros(21,1); v(9) = 1;
assert(abs(H(6,:)*v-n)<1e-12);  % axial rotlet moment
x = [1.1 0.3 0.4;0.5 -0.8 0.9]; normal = x./vecnorm(x,2,2);
angle = 0.43; Q = [cos(angle),-sin(angle),0;sin(angle),cos(angle),0;0,0,1];
[b2,C2] = ring_mixed_basis(Y*Q.',id,Yc*Q.',id,[0 0 0],1,["S","R","T_t","T_n","D","Q_s"]);
G = stokes_mfs_velocity_matrix(basis,x)*C;
G2 = stokes_mfs_velocity_matrix(b2,x*Q.')*C2;
T = stokes_mfs_traction_matrix(basis,x,normal)*C;
T2 = stokes_mfs_traction_matrix(b2,x*Q.',normal*Q.')*C2;
assert(norm(G2-kron(eye(2),Q)*G,'fro')<1e-11);
assert(norm(T2-kron(eye(2),Q)*T,'fro')<1e-10);
% Reflection symmetry tests all supported family-specific component choices.
families = ["S","R","D","T_t","T_n","T_n2","Q_s","S_y"];
for family = families
    [bs,Cs,ws,ks,~,idx] = ring_mixed_basis(Y,id,Yc,id,[0 0 0],1,family,"swirl_free");
    [~,Cf] = ring_mixed_basis(Y,id,Yc,id,[0 0 0],1,family);
    assert(norm(Cs-Cf(:,idx),'fro')<1e-14 && abs(ws.'*ks-1)<1e-14);
    expected = 4;
    if ismember(family,["R","T_n2"]), expected = 3; end
    assert(size(Cs,2)==expected);
    probe = [1.1 0.3 0.4;1.1 -0.3 0.4];
    normal = probe./vecnorm(probe,2,2);
    Gs = stokes_mfs_velocity_matrix(bs,probe)*Cs;
    Ts = stokes_mfs_traction_matrix(bs,probe,normal)*Cs;
    reflection = diag([1,-1,1]);
    assert(norm(Gs(4:6,:)-reflection*Gs(1:3,:),'fro')<1e-11);
    assert(norm(Ts(4:6,:)-reflection*Ts(1:3,:),'fro')<1e-10);
end
% The derivative contributes no force; its source torque is d*(q cross b).
[bs,Cs,~,~,Hs] = ring_mixed_basis(Y,id,Yc,id,[0 0 0],1,"S_y","full");
d = 1-norm(Yc(1,:));
coeff = zeros(size(Cs,2),1); coeff(5) = 1;
lambda = reshape(Cs*coeff,3,[]);
expectedTorque = sum(cross(bs(2).normalDirections,lambda(:,n+1:end).',2),1).'*d;
assert(norm(Hs(1:3,:)*coeff)<1e-14);
assert(norm(Hs(4:6,:)*coeff-expectedTorque)<1e-12 && norm(expectedTorque)>1e-3);
fprintf('ring_mixed_basis self-test passed.\n');
end
