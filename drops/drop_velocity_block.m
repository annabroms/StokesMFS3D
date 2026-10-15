function G = drop_velocity_block(sources, targets, viscosity)
%DROP_VELOCITY_BLOCK Stokeslet velocity block from sources to targets.
% G is 3M-by-3N with interleaved Cartesian components at each node.
% A source strength f produces velocity S*f/viscosity; S already contains
% the factor 1/(8*pi). The exterior viscosity in the drop model is one.

if nargin == 0
    self_test();
    return
end

validateattributes(sources, {'numeric'}, ...
    {'real','finite','2d','ncols',3,'nonempty'}, mfilename, 'sources');
validateattributes(targets, {'numeric'}, ...
    {'real','finite','2d','ncols',3,'nonempty'}, mfilename, 'targets');
validateattributes(viscosity, {'numeric'}, ...
    {'real','finite','scalar','positive'}, mfilename, 'viscosity');

G = stokes_SLP_mat(sources, targets) / viscosity;
assert(all(isfinite(G), 'all'), 'A source coincides with a target.');
end

function self_test()
%SELF_TEST Check the block layout and a radial Stokeslet reference value.
G = drop_velocity_block([0 0 0], [1 0 0; 0 1 0], 2);
assert(isequal(size(G), [6 3]));
reference = diag([1/(8*pi), 1/(16*pi), 1/(16*pi)]);
assert(norm(G(1:3,:) - reference, 'fro') < 1e-14);
reference = diag([1/(16*pi), 1/(8*pi), 1/(16*pi)]);
assert(norm(G(4:6,:) - reference, 'fro') < 1e-14);
fprintf('drop_velocity_block self-test passed.\n');
end
