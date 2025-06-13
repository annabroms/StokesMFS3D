function G = generate_stokes_mat(rin, rout)
%GENERATE_STOKES_MAT Builds dense Stokeslet target-from-source matrix for 3D Stokes flow.
%
%   G = GENERATE_STOKES_MAT(rin, rout)
%
%   Constructs the full dense matrix G representing the (single-layer)
%   Stokeslet kernel in 3D, mapping a vector of point forces at source locations
%   `rin` to velocities at target locations `rout`.
%
%   INPUTS:
%       rin  - N x 3 matrix of source points.
%       rout - M x 3 matrix of target points.
%
%   OUTPUT:
%       G    - 3M x 3N block matrix where each 3x3 block represents the Stokeslet
%              kernel G(x, y) between target x = rout(i,:) and source y = rin(j,:).
%              That is, for each pair (i, j),
%              G(x,y) = (I + r⊗r / |r|^2) / (8π|r|),  where r = x - y.
%
%   NOTES:
%       - This matrix is expensive to construct and store for large N.
%         It is mostly useful for verification, small problems, or dense solvers.
%       - The output format is consistent with applying S * f where f is a
%         stacked 3N x 1 force vector.
%
%   EXAMPLE USAGE:
%       rin = rand(10,3);
%       rout = rand(12,3);
%       S = generate_stokes_mat(rin, rout);
%
%   Anna Broms June, 12, 2025

d1 = bsxfun(@minus,rin(:,1)',rout(:,1)); % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,rin(:,2)',rout(:,2));
d3 = bsxfun(@minus,rin(:,3)',rout(:,3));
r2 = d1.^2+d2.^2+d3.^2;   % dist^2 mat

ir = 1 ./ sqrt(r2);   % 1/r, 
ir3 = ir ./ r2;        	% 1/r^3, 

d1d2 = d1.*d2.*ir3;
d1d3 = d1.*d3.*ir3;
d2d3 = d2.*d3.*ir3;


D1 = zeros(3); D1(1,1) = 1;
D2 = zeros(3); D2(2,2) = 1;
D3 = zeros(3); D3(3,3) = 1;

E1 = zeros(3); E1(1,2) = 1; E1(2,1) = 1;
E2 = zeros(3); E2(1,3) = 1; E2(3,1) = 1;
E3 = zeros(3); E3(3,2) = 1; E3(2,3) = 1;

dxd_ir3 = kron(d1.^2.*ir3,D1) + kron(d2.^2.*ir3,D2) + kron(d3.^2.*ir3,D3) + kron(d1d2,E1) +...
    kron(d1d3,E2) + kron(d2d3,E3);
G = kron(ir,(1/8/pi)*eye(3)) + (1/8/pi)*dxd_ir3; %  incl prefac


end
