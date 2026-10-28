function [strength,diagnostics,fittedVelocity] = solve_wall_image_system( ...
    A,rhs,svdCutoff,label)
%SOLVE_WALL_IMAGE_SYSTEM Direct rectangular resistance solve with size report.
% Uses the same truncated-SVD factors as solve_resistance_wall denseSVD.
% No row or column rescaling, normal equations, or iterative solve is used.
% All rows and columns follow [x1;y1;z1;x2;...]. Relative cutoff only:
% sigma > svdCutoff*sigma_max; absolute cutoff is zero.
% A no-argument call checks a known rank-deficient overdetermined system.

if nargin == 0
    self_test(); strength = []; diagnostics = []; fittedVelocity = []; return
end
[rows,columns] = size(A);
matrixMiB = 8*rows*columns/2^20;
workingMiB = 6*matrixMiB+8*4*columns^2/2^20;
fprintf('%s: %d rows x %d cols, ratio %.3f, A %.2f MiB, SVD budget ~%.2f MiB.\n', ...
    label,rows,columns,rows/columns,matrixMiB,workingMiB);
if rows < 2*columns
    error('wall:ImageOversampling','At least two rows per scalar unknown are required.');
end
if size(rhs,1) ~= rows || size(rhs,2) ~= 1 || mod(rows,3) ~= 0
    error('wall:ImageRhs','RHS must be a column of interleaved 3-vectors.');
end
timer = tic;
[Y,leftVectors,singularValues] = getPseudoFactors(A,svdCutoff,0);
strength = Y*(leftVectors.'*rhs);
solveSeconds = toc(timer);
clear Y leftVectors
fittedVelocity = A*strength;
residual = fittedVelocity-rhs;
diagnostics = struct('rows',rows,'columns',columns,'matrixMiB',matrixMiB, ...
    'workingMiB',workingMiB,'rank',numel(singularValues), ...
    'singularValues',singularValues,'solveSeconds',solveSeconds, ...
    'fittingRelativeTwoNorm',norm(residual)/max(norm(rhs),realmin), ...
    'coefficientNorm',norm(strength),'maxScalarCoefficient',norm(strength,inf));
fprintf('  retained rank %d/%d, fit 2-norm %.3e, |coeff|_2 %.3e, solve %.3fs.\n', ...
    diagnostics.rank,columns,diagnostics.fittingRelativeTwoNorm, ...
    diagnostics.coefficientNorm,solveSeconds);
end

function self_test()
% Repeated columns require a rank-three minimum-norm solution, not rank six.
B = [eye(3);2*eye(3);-eye(3);3*eye(3)];
A = [B,B]; rhs = B*[1;-2;0.5];
[z,d,v] = solve_wall_image_system(A,rhs,1e-12,"direct SVD self-test");
expected = repmat([1;-2;0.5]/2,2,1);
if d.rank ~= 3 || norm(z-expected) > 1e-12 || norm(v-rhs) > 1e-12
    error('wall:ImageSolveTest','Rank-deficient direct solve regression failed.');
end
fprintf('solve_wall_image_system PASS: rank %d, residual %.3e.\n',d.rank,norm(v-rhs));
end
