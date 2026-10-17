function u = direct_image_velocity(basis,coeff,x)
%DIRECT_IMAGE_VELOCITY Sum source-to-target S/R/T/D kernels without dense blocks.
%   u = direct_image_velocity(basis,coeff,x)
%
% basis is the source-body basis from build_image_discretization, coeff{b}
% holds the scaled N_b-by-3 strengths of block b, and x is M-by-3 targets.
% For r=x-y and viscosity one, this applies the same 1/(8*pi) or 1/(4*pi)
% prefactors as stokes_mfs_velocity_matrix. Repeated image-node blocks
% share their M-by-N pair distances. The S family is evaluated by
% getStokesletFlow (direct mode); R, T and D are summed here.
u = zeros(size(x));
vars = struct('fmm',false);
cachedPoints = [];
for b = 1:numel(basis)
    block = basis(b);
    c = coeff{b};
    if block.family == "S"
        flow = getStokesletFlow(reshape(c.',[],1),block.points,x,vars);
        u = u+reshape(flow,3,[]).';
        continue
    end
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
        case "R"
            v = ir3/(8*pi);
            u(:,1) = u(:,1)+sum(v.*(rz.*c(:,2).'-ry.*c(:,3).'),2);
            u(:,2) = u(:,2)+sum(v.*(rx.*c(:,3).'-rz.*c(:,1).'),2);
            u(:,3) = u(:,3)+sum(v.*(ry.*c(:,1).'-rx.*c(:,2).'),2);
        case "T"
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
        otherwise
            error('direct_image_velocity:BadFamily', ...
                'Unknown fundamental-solution family "%s".',block.family);
    end
end
end
