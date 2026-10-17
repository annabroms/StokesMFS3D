function factor = one_body_image_mobility_preconditioner(body,q,cutoff)
%ONE_BODY_IMAGE_MOBILITY_PRECONDITIONER Factor a row-scaled mobility self block.
%   factor = one_body_image_mobility_preconditioner(body,q,cutoff)
%   one_body_image_mobility_preconditioner() checks the block and SVD action.
%
% A = W*(G*(I-C*H)+K*H), where W repeats body.rowScale across components,
% with H mapping source coefficients to total
% force/torque, C=H'/(H*H'), and K mapping rigid motion to boundary speed.
% factor.Y*(factor.U'*mu) applies its truncated pseudoinverse.

if nargin == 0
    self_test();
    factor = struct();
    return
end
if nargin < 3 || isempty(cutoff)
    cutoff = 5e-12;
end
[H,C] = image_mobility_maps(body,q);
G = stokes_mfs_velocity_matrix(body.basis,body.collocation);
K = getKmat(body.collocation,q);
A = G+(K-G*C)*H;
A = repelem(body.rowScale,3).*A;
[Y,U,singularValues] = getPseudoFactors(A,cutoff,0);
assert(~isempty(singularValues),'The SVD cutoff removed every mode.');
factor.H = H;
factor.C = C;
factor.Y = Y;
factor.U = U;
factor.rank = numel(singularValues);
factor.singularValues = singularValues;
factor.sourceCount = size(A,2);
factor.collocationCount = size(body.collocation,1);
end

function self_test()
%SELF_TEST compares the weighted mobility block with its explicit form.
root = fileparts(fileparts(mfilename('fullpath')));
oldPath = path; cleanup = onCleanup(@() path(oldPath));
addpath(fullfile(root,'common','dense'));
body.collocation = [1,0,0;0,1,0;0,0,1;-1,0,0;0,-1,0;0,0,-1];
body.weights = (1:6).';
body.basis = struct('family',"S",'points',0.5*body.collocation, ...
    'normalDirections',zeros(0,3),'scale',1);
for scaled = [true,false]
    if scaled
        body.rowScale = sqrt(body.weights);
    else
        body.rowScale = ones(6,1);
    end
    factor = one_body_image_mobility_preconditioner( ...
        body,[0,0,0],1e-12);
    G = stokes_mfs_velocity_matrix(body.basis,body.collocation);
    K = getKmat(body.collocation,[0,0,0]);
    A = G*(eye(size(G,2))-factor.C*factor.H)+K*factor.H;
    A = repelem(body.rowScale,3).*A;
    rhs = sin((1:size(A,1)).');
    expected = pinv(A,1e-12*norm(A,2))*rhs;
    actual = factor.Y*(factor.U'*rhs);
    assert(norm(expected-actual) < 1e-9*max(1,norm(expected)));
end
fprintf('one_body_image_mobility_preconditioner self-test passed.\n');
end
