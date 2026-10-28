function [lambda_gmres,info] = solve_wall_images_svd(problem,rhs)
%SOLVE_WALL_IMAGES_SVD Solve the global rectangular wall-image system.
% Shared by resistance and mobility when opt.SVDsolve=true. This path uses
% DIRECT complete-wall kernels; no FMM or one-body preconditioner is needed.
% Resistance factors the physical velocity matrix. Mobility factors the
% global recompleted matrix: remove each source body's moment-carrying flow
% at ALL targets and add Kout*Kin' only at that body's own targets.
%
% Rows interleave velocities at all rout; columns interleave coefficients
% in the prepared body/family order. lambda_gmres is the physical source
% for resistance or the auxiliary, unprojected source for mobility, both
% dimensionless. The public solvers perform reconstruction and unit conversion.
% The pseudoinverse keeps sigma>opt.tol*sigmaMax. No identity replacement,
% full force/torque projector, normal equations or rank shifts are used.

nRow = problem.targetOffsets(end); nCol = problem.sourceOffsets(end);
matrixMiB = 8*double(nRow)*double(nCol)/2^20;
fprintf('Global %s SVD: %d x %d, matrix %.1f MiB, estimated SVD work >= %.1f MiB.\n', ...
    problem.mode,nRow,nCol,matrixMiB,6*matrixMiB+32*double(nCol)^2/2^20);
timer = tic;
matrix = zeros(nRow,nCol);
options = struct('method','blake','mu',1,'wall_height',0);
for k = 1:problem.P
    body = problem.body{k};
    cols = problem.sourceOffsets(k)+(1:3*problem.N_total(k));
    Sblock = wall_image_velocity_matrix(body.basis,problem.rout,options);
    if strcmp(problem.mode,'mobility')
        % This applies the main mobility solver's (I-LL) on source columns.
        Sblock = Sblock-(Sblock*body.completion)*body.Kin.';
        rows = problem.targetOffsets(k)+(1:3*problem.M(k));
        Sblock(rows,:) = Sblock(rows,:)+body.Kout*body.Kin.';
    end
    matrix(:,cols) = Sblock;
end
assemblySeconds = toc(timer);
timer = tic;
% UU is stored with singular vectors in rows, as in the one-body solvers.
[Y,UU,Si] = getPseudoFactors(matrix,problem.options.tol,0);
UU = UU.';
lambda_gmres = Y*(UU*rhs);
factorSeconds = toc(timer);
residual = norm(matrix*lambda_gmres-rhs);
if norm(rhs)>0, residual = residual/norm(rhs); end
info = struct('rank',numel(Si),'singularValues',Si,'matrixMiB',matrixMiB, ...
    'systemSize',size(matrix),'linearSystemResidual',residual, ...
    'assemblySeconds',assemblySeconds,'factorizationSeconds',factorSeconds);
end
