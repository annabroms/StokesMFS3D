function [S,check] = wall_images_stresslet(details,checkResolution)
%WALL_IMAGES_STRESSLET Physical force-on-fluid stresslet of wall-image solves.
% Uses TOTAL normalized coefficients in details.lambda, not auxiliary ones.
% The common stresslet unit is mu*a^2, including intrinsic T moments.
% For rigid boundary motion S equals the surface STF traction moment.
% For a residual boundary velocity e, Lorentz reciprocity gives
% S_traction = S - 2*mu*STF(integral_surface e*n' dA).
% This includes complete wall fields: smooth reflected fields have no
% enclosed source moment, but their surface-velocity term need not vanish.
% Optional checkResolution requests this correction on a genuine quadrature
% (not fitting caps). Refine it before interpreting the correction as resolved.
% No arguments runs a dimensional scaling and reciprocal-identity self-test.

if nargin == 0
    self_test(); S = []; check = []; return
end
problem = details.problem;
a = problem.options.radius;
mu = problem.options.viscosity;
grid = struct('centres',problem.q,'body',[problem.body{:}]);
strength = cell(problem.P,1);
for p = 1:problem.P
    cols = problem.sourceOffsets(p)+(1:3*problem.N_total(p));
    strength{p} = details.lambda(cols);
end
S = mu*a^2*apply_KS_transpose(grid,strength);
check = struct('tractionStresslet',[],'correction',[], ...
    'relativeCorrection',NaN,'resolution',[]);
if nargin < 2 || isempty(checkResolution), return; end
surface = setupsurfquad(ellipsoid_param(1,1,1),checkResolution);
correction = zeros(3,3,problem.P);
for p = 1:problem.P
    targets = problem.q(p,:)+surface.x.';
    velocity = reshape(apply_wall_images_velocity( ...
        details.lambda,problem,targets,true),3,[]).';
    U = details.rigidVelocity(6*(p-1)+(1:6));
    U(4:6) = a*U(4:6);
    rigid = reshape(getKmat(targets,problem.q(p,:))*U,3,[]).';
    moment = (velocity-rigid).'*(surface.w(:).*surface.nx.');
    correction(:,:,p) = 2*mu*a^2*stf(moment);
end
check.tractionStresslet = S-correction;
check.correction = correction;
check.relativeCorrection = norm(correction(:))/max(norm(S(:)),realmin);
check.resolution = checkResolution;
end

function value = stf(moment)
% Extract the symmetric trace-free part of a Cartesian moment.
value = (moment+moment.')/2-trace(moment)*eye(3)/3;
end

function self_test()
% Verify physical conversion and the reciprocal identity against analytic traction.
a = 1.7; mu = 2.3; q = [0.2,-0.1,1.4];
block = struct('family',"S",'points',q+[0.12,-0.08,0.04], ...
    'normalDirections',zeros(0,3),'scale',1);
basis = repmat(block,4,1);
families = ["S","T","D","R"];
for b = 1:4
    basis(b).family = families(b);
    if families(b)=="T", basis(b).normalDirections = [0,0,-1]; end
end
lambda = sin((1:12).');
problem = struct('options',struct('radius',a,'viscosity',mu), ...
    'q',q,'P',1,'body',{{struct('basis',basis)}}, ...
    'sourceOffsets',[0;12],'N_total',4);
details = struct('problem',problem,'lambda',lambda);
S = wall_images_stresslet(details);
physicalBasis = basis;
powers = [0,1,2,1]; physicalLambda = lambda;
for b = 1:4
    physicalBasis(b).points = a*basis(b).points;
    cols = 3*(b-1)+(1:3);
    physicalLambda(cols) = mu*a^(1+powers(b))*lambda(cols);
end
physicalGrid = struct('centres',a*q,'body',struct('basis',physicalBasis));
reference = apply_KS_transpose(physicalGrid,{physicalLambda});
assert(norm(S-reference,'fro')<1e-12);

% Arbitrary mixed sources are not rigid: the surface-velocity term is essential.
surface = setupsurfquad(ellipsoid_param(1,1,1),[40,48]);
x = q+surface.x.'; n = surface.nx.';
u = reshape(stokes_mfs_velocity_matrix(basis,x)*lambda,3,[]).';
t = direct_image_traction(basis,lambda,x,n);
grid = struct('centres',q,'body',struct('basis',basis));
source = apply_KS_transpose(grid,{lambda});
traction = stf(surface.x*(surface.w(:).*t));
velocityTerm = 2*stf(u.'*(surface.w(:).*n));
identityError = norm(traction-source+velocityTerm,'fro');
assert(identityError<1e-10);
assert(norm(velocityTerm,'fro')>1e-3);

% A smooth incompressible strain models the regular part of a wall field.
E = diag([1,-0.3,-0.7]);
smoothVelocity = surface.x.'*E.';
smoothTraction = -2*n*E.';
smoothMoment = stf(surface.x*(surface.w(:).*smoothTraction));
smoothCorrection = 2*stf(smoothVelocity.'*(surface.w(:).*n));
assert(norm(smoothMoment+smoothCorrection,'fro')<1e-12);
fprintf('wall_images_stresslet PASS: physical scaling; reciprocal error %.3e.\n', ...
    identityError);
end
