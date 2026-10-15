function coefficients = applyOneBodyPinv(boundaryData, factors)
%APPLYONEBODYPINV Apply the retained local SVD in stable factor order.
% First take U_r'*boundaryData, then divide by sigma_r, then apply V_r.
% Column scales are undone last, giving physical [alpha;gamma].

if nargin == 0
    self_test();
    return
end

assert(iscolumn(boundaryData) && ...
    numel(boundaryData) == size(factors.leftTranspose,2), ...
    'Boundary data must have one entry per local row.');
leftCoordinates = factors.leftTranspose*boundaryData;
singularCoordinates = leftCoordinates./factors.singularValues;
scaledCoefficients = factors.rightVectors*singularCoordinates;
coefficients = scaledCoefficients./factors.columnScales;
end

function self_test()
% Check a tall reference solve and the order of optional column unscaling.
matrix = [2 0;0 1;0 0];
[left,singular,right] = svd(matrix,'econ');
factors.leftTranspose = left.';
factors.singularValues = diag(singular);
factors.rightVectors = right;
factors.columnScales = ones(2,1);
eta = [2;2;7];
assert(norm(applyOneBodyPinv(eta,factors)-[1;2]) < 1e-14);

factors.columnScales = [2;0.5];
[left,singular,right] = svd(matrix./factors.columnScales.','econ');
factors.leftTranspose = left.';
factors.singularValues = diag(singular);
factors.rightVectors = right;
assert(norm(applyOneBodyPinv(eta,factors)-[1;2]) < 1e-14);
fprintf('applyOneBodyPinv self-test passed.\n');
end
