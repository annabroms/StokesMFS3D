function K = getKmat(r, q)
%K = GETKMAT(r, q) returns a matrix K that maps a vector of rigid body
%   velocities (translations and rotations) into the corresponding
%   pointwise velocities at a set of surface points. The transpose, K',
%   maps a force density on surfaces to the net force and torque on the
%   bodies. To be used in MFS or the rigid multiblob method.
%
%   INPUTS:
%     r - (NP x 3) array of 3D coordinates of points on the surfaces of P particles.
%         These are assumed to be ordered consecutively per particle, with N the 
%         number of points per particle (assumed equal for everone).
%
%     q - (P x 3) array of 3D coordinates for the centers of the P spherical particles.
%
%   NOTE: Most commonly, this function is used for one particle at the time.
%
%   OUTPUT:
%     K - (3N x 6P) matrix such that K * U maps rigid body velocities
%         to surface point velocities:
%           - U is a 6P vector [u1; omega1; u2; omega2; ...], where each
%             ui (3x1) is the translational velocity and omegai (3x1) is the angular
%             velocity of particle i.
%           - K * U yields the stacked velocities (3N x 1) of all surface points.
%
%   INTERNALS:
%     For each particle i:
%       - A 3x3 identity block in J corresponds to translation (u_i).
%       - The skew-symmetric matrix (via cross_mat) maps angular velocity 
%         to pointwise rotational contribution: omega x (r - q_i).
%       - Each 3N x 6 block of K is [I; cross_mat(r - q_i)], 
%
%   REFERENCE:
%     The sign convention for the cross product matrix follows the right-hand rule,
%     as per: https://en.wikipedia.org/wiki/Cross_product

% Construct helper function to get cross product matrix (used for omega × x)

assert(size(r,2)==3,"r must be transposed")
assert(size(q,2)==3,"q must be transposed")

cross_mat = @(x) -[0    -x(3)  x(2);
                   x(3)  0    -x(1);
                  -x(2)  x(1)   0];

NP = size(r,1);       % Total number of surface points
P = size(q,1);       % Number of particles
N = NP / P;          % Number of points per particle (assumed equal)

% Initialize K
K = zeros(3*NP, 6*P); % Final matrix to fill

% Create repeated identity blocks for translation
J = repmat(eye(3), N, 1);  % 3*N x 3

for i = 1:P
    B = zeros(3*N, 3);     % Placeholder for rotational block
    
    % Compute (r - q_i) for each point on particle i
    for k = 1:N
        idx = (i-1)*N + k;         % global index for point
        B(3*(k-1)+1:3*k,:) = cross_mat(r(idx,:) - q(i,:));
    end

    % Assemble 3N x 6 block: [I; skew]
    K_block = [J, B];

    % Place block into K
    row_range = (i-1)*3*N+1 : i*3*N;
    col_range = (i-1)*6+1   : i*6;
    K(row_range, col_range) = K_block;
end

end



