function [H,C] = image_mobility_maps(body,q)
%IMAGE_MOBILITY_MAPS Map all source coefficients to total force and torque.
%   [H,C] = image_mobility_maps(body,q)
%   image_mobility_maps() checks mixed Stokeslet/rotlet moments.
%
% H is 6-by-n and C = H'/(H*H') is n-by-6. S_y columns have zero force
% and source torque d*(q cross b). A coefficient vector z has
% force/torque H*z; C*F is its minimum-norm completion. The total-moment
% projector is L*z = C*(H*z). Blocks follow body.basis order and each block
% has interleaved [x1;y1;z1;x2;...] components. Viscosity is one.

if nargin == 0
    self_test();
    H = []; C = [];
    return
end
n = 3*sum(arrayfun(@(b) size(b.points,1),body.basis));
H = zeros(6,n);
offset = 0;
for b = 1:numel(body.basis)
    block = body.basis(b);
    count = size(block.points,1);
    cols = offset+(1:3*count);
    if block.family == "S"
        H(:,cols) = block.scale*getKmat(block.points,q).';
    elseif block.family == "R"
        H(4:6,cols) = block.scale*repmat(eye(3),1,count);
    elseif block.family == "S_y"
        assert(isequal(size(block.normalDirections),[count 3]) && ...
            isequal(size(block.derivativeDistances),[count 1]));
        for j = 1:count
            qj = block.normalDirections(j,:);
            crossMatrix = [0,-qj(3),qj(2);qj(3),0,-qj(1); ...
                -qj(2),qj(1),0];
            H(4:6,offset+3*j-2:offset+3*j) = ...
                block.scale*block.derivativeDistances(j)*crossMatrix;
        end
    end
    offset = offset+3*count;
end
C = H'/(H*H');
end

function self_test()
%SELF_TEST checks scaled moments, total cancellation, and completion.
q = [0,0,0];
base = struct('family',"S",'points', ...
    [1,0,0;-1,0,0;0,1,0;0,0,1], ...
    'normalDirections',zeros(0,3),'scale',1);
rotlet = base;
rotlet.family = "R";
rotlet.points = [0.3,0,0];
rotlet.scale = 2;
dipole = base;
dipole.family = "D";
dipole.points = [0.2,0,0];
stresslet = dipole;
stresslet.family = "T";
body.basis = [base;rotlet;dipole;stresslet];
quadrupole = dipole;
quadrupole.family = "Q_s";
quadrupole.normalDirections = [1 0 0];
body.basis(end+1) = quadrupole;
[H,C] = image_mobility_maps(body,q);
assert(isequal(size(H),[6,24]) && isequal(size(C),[24,6]));
assert(all(H(:,16:24) == 0,'all'));
z = zeros(24,1);
z(2) = 1; z(5) = -1; z(15) = -1; % R has scale two.
assert(norm(H*z) < 1e-14);
f = [1;2;3;4;5;6];
assert(norm(H*(C*f)-f) < 1e-13);
L = C*H;
assert(norm(L*L-L,'fro') < 1e-12);
assert(norm(H*(z-L*z)) < 1e-13);
base.derivativeDistances = zeros(0,1);
derivative = base;
derivative.family = "S_y";
derivative.points = [0.2 0 0];
derivative.normalDirections = [0 1 0];
derivative.derivativeDistances = 0.3;
body.basis = [base(1);derivative];
[Hd,~] = image_mobility_maps(body,q);
probe = zeros(size(Hd,2),1);
probe(end) = 2;
assert(norm(Hd(1:3,:)*probe)<1e-14 && ...
    norm(Hd(4:6,:)*probe-[0.6;0;0])<1e-14);
fprintf('image_mobility_maps self-test passed.\n');
end
