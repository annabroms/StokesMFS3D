function u = direct_image_velocity(basis,coeff,x)
%DIRECT_IMAGE_VELOCITY Sum source-to-target image kernels without dense blocks.
%   u = direct_image_velocity(basis,coeff,x)
%
% basis is the source-body basis from build_image_discretization, coeff{b}
% holds the scaled N_b-by-3 strengths of block b, and x is M-by-3 targets.
% For r=x-y and viscosity one, this applies the same 1/(8*pi) or 1/(4*pi)
% prefactors as stokes_mfs_velocity_matrix. Repeated image-node blocks
% share their M-by-N pair distances; all families are summed directly.
if nargin == 0
    self_test(); u = []; return
end
u = zeros(size(x));
cachedPoints = [];
for b = 1:numel(basis)
    block = basis(b);
    c = coeff{b};
    if ~isequal(block.points,cachedPoints)
        y = block.points;
        cachedPoints = y;
        rx = x(:,1)-y(:,1).';
        ry = x(:,2)-y(:,2).';
        rz = x(:,3)-y(:,3).';
        r2 = rx.^2+ry.^2+rz.^2;
        ir = 1./sqrt(r2);
        ir3 = ir./r2;
    end
    switch block.family
        case "S"
            rc = rx.*c(:,1).'+ry.*c(:,2).'+rz.*c(:,3).';
            u(:,1) = u(:,1)+(ir*c(:,1)+sum(rx.*rc.*ir3,2))/(8*pi);
            u(:,2) = u(:,2)+(ir*c(:,2)+sum(ry.*rc.*ir3,2))/(8*pi);
            u(:,3) = u(:,3)+(ir*c(:,3)+sum(rz.*rc.*ir3,2))/(8*pi);
        case "R"
            v = ir3/(8*pi);
            u(:,1) = u(:,1)+sum(v.*(rz.*c(:,2).'-ry.*c(:,3).'),2);
            u(:,2) = u(:,2)+sum(v.*(rx.*c(:,3).'-rz.*c(:,1).'),2);
            u(:,3) = u(:,3)+sum(v.*(ry.*c(:,1).'-rx.*c(:,2).'),2);
        case {"T","T_n","T_n2","T_t"}
            n = block.normalDirections;
            assert(isequal(size(n),size(c)));
            rn = rx.*n(:,1).'+ry.*n(:,2).'+rz.*n(:,3).';
            rc = rx.*c(:,1).'+ry.*c(:,2).'+rz.*c(:,3).';
            v = (-6/(8*pi))*(rn.*rc.*ir3./r2);
            u(:,1) = u(:,1)+sum(rx.*v,2);
            u(:,2) = u(:,2)+sum(ry.*v,2);
            u(:,3) = u(:,3)+sum(rz.*v,2);
        case "D"
            rc = rx.*c(:,1).'+ry.*c(:,2).'+rz.*c(:,3).';
            v = (3/(4*pi))*(ir3./r2).*rc;
            u(:,1) = u(:,1)-(ir3*c(:,1))/(4*pi)+sum(rx.*v,2);
            u(:,2) = u(:,2)-(ir3*c(:,2))/(4*pi)+sum(ry.*v,2);
            u(:,3) = u(:,3)-(ir3*c(:,3))/(4*pi)+sum(rz.*v,2);
        case "Q_s"
            u = u+stokes_axial_quadrupole_field( ...
                block.points,x,block.normalDirections,c,"velocity");
        otherwise
            error('direct_image_velocity:BadFamily', ...
                'Unknown fundamental-solution family "%s".',block.family);
    end
end
end

function self_test()
% Compare all direct family sums with the mixed dense velocity block.
b = struct('family',"S",'points',[0 0 0], ...
    'normalDirections',zeros(0,3),'scale',1);
basis = repmat(b,8,1);
families = ["S","R","T","D","T_n","T_t","T_n2","Q_s"];
for k = 1:numel(families)
    basis(k).family = string(families(k));
    basis(k).points = [0.1*k,0.03*k,-0.02*k];
    basis(k).scale = 0.6+0.2*k;
    if ismember(k,[3 5 6 7 8]), basis(k).normalDirections = [0 1 0]; end
end
z = sin((1:3*numel(families)).');
x = [1 0.2 0.3;0.3 1 -0.2];
c = image_source_coefficients(basis,z);
actual = direct_image_velocity(basis,c,x);
expected = reshape(stokes_mfs_velocity_matrix(basis,x)*z,3,[]).';
assert(norm(actual-expected,inf)<1e-12);
fprintf('direct_image_velocity self-test passed.\n');
end
