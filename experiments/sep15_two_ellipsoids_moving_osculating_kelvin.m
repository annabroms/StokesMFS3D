function results = sep15_two_ellipsoids_moving_osculating_kelvin(visualise,options)
%SEP15_TWO_ELLIPSOIDS_MOVING_OSCULATING_KELVIN Moving sphere image clouds.
%
%   results = sep15_two_ellipsoids_moving_osculating_kelvin;
%   results = sep15_two_ellipsoids_moving_osculating_kelvin(false);
%   results = sep15_two_ellipsoids_moving_osculating_kelvin(true, ...
%       struct('nTheta',96,'nRadial',16,'maxCycles',250));
%   checks = sep15_two_ellipsoids_moving_osculating_kelvin('self_test');
%
%   Two close, rotated, triaxial ellipsoids are reflected into alternately.
%   Start with a deterministic sample of the FILLED focal ellipse of body 2:
%       x^2/(a^2-c^2) + y^2/(b^2-c^2) <= 1, z = 0
%   in its body frame, where a>b>c. This is a two-dimensional disk, NOT a
%   volume or the confocal ellipsoid at a nondegenerate parameter value.
%   The existing discretize_confocal_ellipse utility supplies uniform-area
%   interior rings; a centre and boundary nodes complete the plotting mesh.
%   A finite discretisation samples, rather than represents, all points.
%
%   Every query is projected onto the other ellipsoid afresh. At that
%   closest point p, compute outward normal n and positive principal
%   curvatures kmin,kmax. A tangent sphere cannot match both curvatures
%   unless p is umbilic. Compare three DISTINCT, deterministic models:
%       mean: R = 2/(kmin+kmax)    (mean curvature, not mean radius)
%       min:  R = 1/kmin           (largest principal radius)
%       max:  R = 1/kmax           (smallest principal radius).
%   The exterior displacement q-p is normal, so it does not select a tangent
%   direction. Retaining BOTH curvatures instead requires an anisotropic
%   local model, e.g. the osculating paraboloid in principal coordinates
%       z = -(kmin*u^2+kmax*v^2)/2.
%   Such a paraboloid does not have a standard sphere Kelvin inversion.
%   Its centre is C=p-R*n and Kelvin inversion is
%       I(q)=C+R^2*(q-C)/norm(q-C)^2.
%   As q=p+d*n for exterior closest-point projection, evaluate this stably
%   as I(q)=p-[R*d/(R+d)]*n. Radii and contacts MOVE with each query.
%   This is a geometric surrogate, not an exact ellipsoid Schwarz/Kelvin
%   reflection or a Stokes image system. A tangent sphere need not be
%   contained in the ellipsoid; image-domain violations are reported.
%   There are no MFS strengths, collocation operators, or PDE solves here.
%
%   One cycle consists of I1 followed by I2. Save body-1 and body-2 images
%   separately: convergence is a two-cycle, not a point common to both.
%   For each model, plot the advected disk and ring curves, cloud diameter,
%   triangular area, three principal RMS widths, and endpoint residuals.
%   These are images of the initial parameterised mesh, not contour lines
%   of a PDE solution. Connectivity is retained without remeshing; triangles
%   interpolate mapped nodes, NOT an exact map of triangle interiors.
%   A curved 2D sheet can have THREE nonzero global PCA widths; PCA rank is
%   not an intrinsic-dimension test. Collapse is tolerance-dependent and
%   is not an exact finite-iteration loss of dimension.
%
%   GEOMETRIC INVARIANT: fixed interior image pairs lie on the closest-gap
%   normal. At a two-cycle both points and their closest contacts are on
%   the same line; the contacts have opposing outward normals. If images
%   are inside their respective disjoint strictly convex bodies, these
%   contacts are the unique minimum-distance pair. A fixed local radius
%   rule therefore has at most one interior fixed pair. This does not prove
%   that every starting point converges, or rule out escaping/non-fixed
%   dynamics. Compare with the analytic sphere pair at those contacts and
%   the numerical Jacobian of the FULL moving-sphere two-step map.
%
%   A separate directional-radius sweep at the gap predicts axial intervals
%   as the normal-section direction varies. It is a FAMILY of models, not
%   evidence of start dependence of any one of the three maps above.
%
%   Options (all optional): nTheta=48, nRadial=8, gap=0.08,
%       maxCycles=180, orbitTol=1e-12, endpointTol=1e-7,
%       collapseTols=[1e-3 1e-6 1e-9],
%       cyclesToShow=[0 1 3 8 16 32 64 120], directionalCount=361.
%   Increase sampling and tighten tolerances before interpreting a limit
%   cloud. 'self_test' runs small geometry/algebra tests WITHOUT figures
%   and does not run the cloud experiment.
%
%   References: D. Eberly, DistancePointEllipseEllipsoid.pdf,
%   https://www.geometrictools.com/Documentation/DistancePointEllipseEllipsoid.pdf
%   Principal curvatures and normal sections:
%   https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node30.html

if nargin >= 1 && (ischar(visualise) || isstring(visualise))
    if strcmp(visualise,'self_test')
        results = self_test;
        return;
    end
    error('moving_kelvin:InvalidMode','The only text mode is self_test.');
end
if nargin < 1 || isempty(visualise)
    visualise = true;
end
if nargin < 2 || isempty(options)
    options = struct;
end
validateattributes(visualise,{'logical','numeric'},{'scalar','binary'});
defaults = struct('nTheta',48,'nRadial',8,'gap',0.08, ...
    'maxCycles',180,'orbitTol',1e-12,'endpointTol',1e-7, ...
    'collapseTols',[1e-3 1e-6 1e-9], ...
    'cyclesToShow',[0 1 3 8 16 32 64 120],'directionalCount',361);
assert(isstruct(options) && isscalar(options),'options must be a scalar struct.');
names = fieldnames(options);
for k = 1:numel(names)
    assert(isfield(defaults,names{k}),'Unknown option: %s.',names{k});
    defaults.(names{k}) = options.(names{k});
end
options = defaults;
integerFields = {'nTheta','nRadial','maxCycles','directionalCount'};
for k = 1:numel(integerFields)
    validateattributes(options.(integerFields{k}),{'numeric'}, ...
        {'real','finite','scalar','positive','integer'});
end
assert(options.nTheta >= 8 && options.directionalCount >= 3, ...
    'Use at least 8 angular nodes and 3 directional samples.');
positiveFields = {'gap','orbitTol','endpointTol'};
for k = 1:numel(positiveFields)
    validateattributes(options.(positiveFields{k}),{'numeric'}, ...
        {'real','finite','scalar','positive'});
end
validateattributes(options.collapseTols,{'numeric'}, ...
    {'real','finite','vector','nonempty','positive'});
validateattributes(options.cyclesToShow,{'numeric'}, ...
    {'real','finite','vector','nonempty','integer','nonnegative'});

% Invariant: 3-by-N point columns in world coordinates; normals are outward.
% Separate support planes certify the exact gap without optimisation.
[bodies,gapGeometry] = make_geometry(options.gap);
[seeds,mesh] = make_focal_seeds(bodies(2),options);
rules = {'mean','min','max'};
runs = cell(1,numel(rules));
for ruleIndex = 1:numel(rules)
    reference = gap_reference(bodies,gapGeometry,rules{ruleIndex});
    runs{ruleIndex} = iterate_cloud(seeds,mesh,bodies, ...
        rules{ruleIndex},reference,options);
end
directional = directional_gap_family(bodies,gapGeometry,options.directionalCount);
results = struct('options',options,'bodies',bodies,'gap',gapGeometry, ...
    'seeds',seeds,'mesh',mesh,'rules',{rules},'runs',{runs}, ...
    'directional_family',directional);
report_results(results);
if visualise
    plot_results(results);
end
end

function [bodies,gapGeometry] = make_geometry(gap)
% Set exact facing support points for two asymmetric low-curvature patches.
axes1 = [1.70;1.25;0.72];
axes2 = [1.55;1.15;0.76];
t = 6*pi/180;
s = -4*pi/180;
% Active, right-handed rotations; Q columns are world-space body axes.
Ry = [cos(t),0,sin(t);0,1,0;-sin(t),0,cos(t)];
Rx = [1,0,0;0,cos(s),-sin(s);0,sin(s),cos(s)];
t = -9*pi/180;
s = 5*pi/180;
Rz = [cos(t),-sin(t),0;sin(t),cos(t),0;0,0,1];
Ry2 = [cos(s),0,sin(s);0,1,0;-sin(s),0,cos(s)];
bodies(1) = struct('axes',axes1,'Q',Ry*Rx,'c',zeros(3,1));
bodies(2) = struct('axes',axes2,'Q',Rz*Ry2,'c',zeros(3,1));
n = [0;0;1];
p1 = -gap/2*n;
p2 = gap/2*n;
for k = 1:2
    direction = (3-2*k)*n;
    localDirection = bodies(k).Q.'*direction;
    support = bodies(k).Q*(bodies(k).axes.^2.*localDirection) / ...
        norm(bodies(k).axes.*localDirection);
    desiredPoint = p1;
    if k == 2
        desiredPoint = p2;
    end
    bodies(k).c = desiredPoint-support;
end
gapGeometry = struct('distance',gap,'normal',n, ...
    'contacts',[p1,p2],'midpoint',(p1+p2)/2);
end

function [seeds,mesh] = make_focal_seeds(body,options)
% Reuse the filled-focal-ellipse nodes and keep their ring connectivity.
a = body.axes(1);
b = body.axes(2);
c = body.axes(3);
focalAxes = sqrt([a^2-c^2;b^2-c^2]);
[interior,r,theta] = discretize_confocal_ellipse( ...
    a,b,c,options.nTheta,options.nRadial,2);
boundaryTheta = (0:options.nTheta-1)'*2*pi/options.nTheta;
boundary = [focalAxes(1)*cos(boundaryTheta), ...
    focalAxes(2)*sin(boundaryTheta),zeros(options.nTheta,1)];
localPoints = [zeros(1,3);interior;boundary];
normalised = localPoints(:,1:2)./focalAxes.';
triangles = delaunay(normalised(:,1),normalised(:,2));
seeds = body.c+body.Q*localPoints.';
assert(all(ellipsoid_level(seeds,body) < 1), ...
    'The entire focal disk must be strictly inside its ellipsoid.');
ringIndices = cell(1,options.nRadial+1);
for k = 1:options.nRadial
    ringIndices{k} = 1+(k-1)*options.nTheta+(1:options.nTheta);
end
ringIndices{end} = 1+size(interior,1)+(1:options.nTheta);
mesh = struct('triangles',triangles,'local_points',localPoints.', ...
    'normalised_points',normalised.','focal_axes',focalAxes, ...
    'radial_coordinate',[0;r;ones(options.nTheta,1)], ...
    'azimuth',[0;theta;boundaryTheta],'ring_indices',{ringIndices});
end

function level = ellipsoid_level(points,body)
% Test actual ellipsoid membership, rather than tangent-sphere membership.
local = body.Q.'*(points-body.c);
level = sum((local./body.axes).^2,1);
end

function [contact,normal,curvatures,projectionError] = closest_jet(query,body)
% Project exterior queries via a monotone scalar root and compute both curvatures.
localQuery = body.Q.'*(query-body.c);
assert(size(query,1) == 3 && all(isfinite(query(:))), ...
    'Expected finite real 3-by-N query columns.');
assert(all(sum((localQuery./body.axes).^2,1) > 1), ...
    ['A query entered the ellipsoid being reflected into. Interior closest ' ...
     'projection can be nonunique and is intentionally not continued.']);
axisSquared = body.axes.^2;
lower = zeros(1,size(query,2));
upper = sqrt(sum((body.axes.*localQuery).^2,1));
% f(lambda)=sum((a_i*q_i/(a_i^2+lambda))^2)-1 decreases on [0,infinity).
% This upper bound has f(upper)<=0, so bisection is globally bracketed.
for iteration = 1:80
    lambda = (lower+upper)/2;
    value = sum((body.axes.*localQuery./(axisSquared+lambda)).^2,1)-1;
    positive = value > 0;
    lower(positive) = lambda(positive);
    upper(~positive) = lambda(~positive);
end
lambda = (lower+upper)/2;
localContact = axisSquared.*localQuery./(axisSquared+lambda);
contact = body.c+body.Q*localContact;
inverseSquared = 1./axisSquared;
halfGradient = inverseSquared.*localContact;
gradientNorm = sqrt(sum(halfGradient.^2,1));
localNormal = halfGradient./gradientNorm;
normal = body.Q*localNormal;

% The positive shape operator is T'*A*T / norm(A*(p-c)), A=Q*diag(1/a_i^2)*Q'.
% Trace and determinant of its 2-by-2 tangent restriction avoid choosing
% principal directions, which are indeterminate at umbilics.
traceK = (sum(inverseSquared)- ...
    sum(inverseSquared.*localNormal.^2,1))./gradientNorm;
cofactors = [inverseSquared(2)*inverseSquared(3); ...
    inverseSquared(1)*inverseSquared(3); ...
    inverseSquared(1)*inverseSquared(2)];
detK = sum(cofactors.*localNormal.^2,1)./gradientNorm.^2;
discriminant = traceK.^2-4*detK;
assert(all(discriminant >= -200*eps*traceK.^2), ...
    'Invalid principal-curvature discriminant.');
% Only negative roundoff in the exact nonnegative discriminant is clipped.
kmax = (traceK+sqrt(max(discriminant,0)))/2;
kmin = detK./kmax;
curvatures = [kmin;kmax];
displacement = query-contact;
distance = sum(displacement.*normal,1);
tangentError = displacement-normal.*distance;
projectionError = [max(abs(sum((localContact./body.axes).^2,1)-1)); ...
    max(sqrt(sum(tangentError.^2,1))./max(1,sqrt(sum(displacement.^2,1))))];
assert(all(projectionError < 2e-12), ...
    'Closest-point boundary or normality check failed.');
end

function radius = choose_radius(curvatures,rule)
% Keep the sphere-selection convention explicit and deterministic.
switch rule
    case 'mean'
        radius = 2./sum(curvatures,1);
    case 'min'
        radius = 1./curvatures(1,:);
    case 'max'
        radius = 1./curvatures(2,:);
    otherwise
        error('moving_kelvin:UnknownRule','Unknown curvature rule %s.',rule);
end
end

function [image,info] = moving_step(query,body,rule)
% Perform one Kelvin inversion in the newly selected tangent sphere.
[contact,normal,curvatures,projectionError] = closest_jet(query,body);
radius = choose_radius(curvatures,rule);
distance = sum((query-contact).*normal,1);
assert(all(distance > 0),'Exterior closest-point distances must be positive.');
depth = radius.*distance./(radius+distance);
image = contact-normal.*depth;
info = struct('contact',contact,'normal',normal,'curvatures',curvatures, ...
    'radius',radius,'centre',contact-normal.*radius,'distance',distance, ...
    'image_level',ellipsoid_level(image,body),'projection_error',projectionError);
end

function curvatures = curvature_at_contact(contact,body)
% Evaluate the tangent shape operator directly for reference points and tests.
A = body.Q*diag(1./body.axes.^2)*body.Q.';
gradient = A*(contact-body.c);
n = gradient/norm(gradient);
[~,index] = min(abs(n));
coordinate = zeros(3,1);
coordinate(index) = 1;
t1 = coordinate-n*dot(n,coordinate);
t1 = t1/norm(t1);
T = [t1,cross(n,t1)];
K = T.'*A*T/norm(gradient);
curvatures = sort(eig(K));
end

function pair = sphere_gap_pair(radius1,radius2,gapGeometry)
% Calculate the attracting pair of fixed-gap tangent spheres on their axis.
gap = gapGeometry.distance;
discriminant = gap.*(2*radius1+gap).*(2*radius2+gap).* ...
    (2*radius1+2*radius2+gap);
offset = -gap.*(radius1-radius2);
denominator = 2*(radius1+radius2+gap);
z1 = (offset-sqrt(discriminant))./denominator;
z2 = (offset+sqrt(discriminant))./denominator;
pair = struct('body1',gapGeometry.midpoint+gapGeometry.normal*z1, ...
    'body2',gapGeometry.midpoint+gapGeometry.normal*z2, ...
    'radii',[radius1;radius2]);
end

function reference = gap_reference(bodies,gapGeometry,rule)
% Verify the analytic gap pair against the full state-dependent maps.
k1 = curvature_at_contact(gapGeometry.contacts(:,1),bodies(1));
k2 = curvature_at_contact(gapGeometry.contacts(:,2),bodies(2));
reference = sphere_gap_pair(choose_radius(k1,rule), ...
    choose_radius(k2,rule),gapGeometry);
mapped1 = moving_step(reference.body2,bodies(1),rule);
mapped2 = moving_step(reference.body1,bodies(2),rule);
reference.residual = max(norm(mapped1-reference.body1), ...
    norm(mapped2-reference.body2));
reference.inside = [ellipsoid_level(reference.body1,bodies(1)); ...
    ellipsoid_level(reference.body2,bodies(2))] < 1;
reference.curvatures = [k1,k2];
reference.linearisation = moving_jacobian(reference.body2,bodies,rule);
end

function data = moving_jacobian(point,bodies,rule)
% Differentiate the FULL moving-sphere composition, not frozen-sphere inversion.
step = 1e-5*max(1,norm(point));
jacobians = zeros(3,3,2);
for refinement = 1:2
    h = step/2^(refinement-1);
    for direction = 1:3
        offset = zeros(3,1);
        offset(direction) = h;
        plus1 = moving_step(point+offset,bodies(1),rule);
        minus1 = moving_step(point-offset,bodies(1),rule);
        plus2 = moving_step(plus1,bodies(2),rule);
        minus2 = moving_step(minus1,bodies(2),rule);
        jacobians(:,direction,refinement) = (plus2-minus2)/(2*h);
    end
end
J = jacobians(:,:,2);
eigenvalues = eig(J);
data = struct('jacobian',J,'eigenvalues',eigenvalues, ...
    'singular_values',svd(J),'spectral_radius',max(abs(eigenvalues)), ...
    'step',step/2,'step_refinement_difference', ...
    norm(jacobians(:,:,1)-J,'fro'));
end

function run = iterate_cloud(seeds,mesh,bodies,rule,reference,options)
% Advance a common seed mesh and diagnose collapse separately for each parity.
numberOfSeeds = size(seeds,2);
body1History = nan(3,numberOfSeeds,options.maxCycles+1);
body2History = nan(3,numberOfSeeds,options.maxCycles+1);
body2History(:,:,1) = seeds;
diameters = nan(options.maxCycles+1,2);
areas = nan(options.maxCycles+1,2);
rmsWidths = nan(options.maxCycles+1,3,2);
[diameters(1,2),rmsWidths(1,:,2),areas(1,2)] = cloud_metrics(seeds,mesh.triangles);
cycleChanges = nan(options.maxCycles,1);
firstConvergedCycle = nan(1,numberOfSeeds);
outsideCount = zeros(options.maxCycles,2);
maximumLevels = nan(options.maxCycles,2);
representativeCurvatures = nan(2,2,options.maxCycles);
representativeRadii = nan(2,options.maxCycles);
representativeContacts = nan(3,2,options.maxCycles);
projectionErrors = zeros(2,1);
body2Images = seeds;
body1Images = nan(size(seeds));
for cycle = 1:options.maxCycles
    [new1,info1] = moving_step(body2Images,bodies(1),rule);
    [new2,info2] = moving_step(new1,bodies(2),rule);
    change = sqrt(sum((new2-body2Images).^2,1));
    if cycle > 1
        change = max(change,sqrt(sum((new1-body1Images).^2,1)));
    end
    cycleChanges(cycle) = max(change);
    newlyConverged = cycle > 1 & isnan(firstConvergedCycle) & change <= options.orbitTol;
    firstConvergedCycle(newlyConverged) = cycle;
    body1Images = new1;
    body2Images = new2;
    body1History(:,:,cycle+1) = new1;
    body2History(:,:,cycle+1) = new2;
    [diameters(cycle+1,1),rmsWidths(cycle+1,:,1),areas(cycle+1,1)] = ...
        cloud_metrics(new1,mesh.triangles);
    [diameters(cycle+1,2),rmsWidths(cycle+1,:,2),areas(cycle+1,2)] = ...
        cloud_metrics(new2,mesh.triangles);
    outsideCount(cycle,:) = [nnz(info1.image_level > 1+2e-12), ...
        nnz(info2.image_level > 1+2e-12)];
    maximumLevels(cycle,:) = [max(info1.image_level),max(info2.image_level)];
    representativeCurvatures(:,:,cycle) = ...
        [info1.curvatures(:,end),info2.curvatures(:,end)];
    representativeRadii(:,cycle) = [info1.radius(end);info2.radius(end)];
    representativeContacts(:,:,cycle) = [info1.contact(:,end),info2.contact(:,end)];
    projectionErrors = max(projectionErrors, ...
        max(info1.projection_error,info2.projection_error));
    if cycle >= 2 && cycleChanges(cycle) <= options.orbitTol
        break;
    end
end
lastCycle = cycle;
mapped1 = moving_step(body2Images,bodies(1),rule);
mapped2 = moving_step(body1Images,bodies(2),rule);
pairResiduals = max(sqrt(sum((mapped1-body1Images).^2,1)), ...
    sqrt(sum((mapped2-body2Images).^2,1)));
joint = [body1Images;body2Images];
[clusterLabels,clusterRepresentatives] = endpoint_clusters(joint,options.endpointTol);
collapseCycles = nan(numel(options.collapseTols),2);
for parity = 1:2
    for toleranceIndex = 1:numel(options.collapseTols)
        threshold = options.collapseTols(toleranceIndex);
        first = find(diameters(1:lastCycle+1,parity) <= threshold,1);
        if ~isempty(first)
            collapseCycles(toleranceIndex,parity) = first-1;
        end
    end
end
referenceErrors = [max(sqrt(sum((body1Images-reference.body1).^2,1))); ...
    max(sqrt(sum((body2Images-reference.body2).^2,1)))];
converged = cycleChanges(lastCycle) <= options.orbitTol && ...
    all(pairResiduals <= 10*options.orbitTol);
pointCollapse = converged && all(diameters(lastCycle+1,:) <= options.endpointTol);
fixedSetStatus = 'Unresolved: no intrinsic dimension is inferred from PCA.';
if pointCollapse
    fixedSetStatus = 'One point per parity, resolved only to the requested tolerances.';
end
run = struct('rule',rule,'cycles',0:lastCycle, ...
    'body1_images',body1History(:,:,1:lastCycle+1), ...
    'body2_images',body2History(:,:,1:lastCycle+1), ...
    'diameters',diameters(1:lastCycle+1,:), ...
    'rms_widths',rmsWidths(1:lastCycle+1,:,:), ...
    'triangular_areas',areas(1:lastCycle+1,:), ...
    'cycle_changes',cycleChanges(1:lastCycle), ...
    'first_converged_cycle',firstConvergedCycle,'pair_residuals',pairResiduals, ...
    'outside_counts',outsideCount(1:lastCycle,:), ...
    'maximum_ellipsoid_levels',maximumLevels(1:lastCycle,:), ...
    'representative_curvatures',representativeCurvatures(:,:,1:lastCycle), ...
    'representative_radii',representativeRadii(:,1:lastCycle), ...
    'representative_contacts',representativeContacts(:,:,1:lastCycle), ...
    'projection_errors',projectionErrors, ...
    'cluster_labels',clusterLabels,'cluster_representatives',clusterRepresentatives, ...
    'collapse_cycles',collapseCycles,'reference',reference, ...
    'reference_errors',referenceErrors,'converged',converged, ...
    'point_collapse_at_tolerance',pointCollapse,'fixed_set_status',fixedSetStatus);
end

function [diameter,widths,area] = cloud_metrics(points,triangles)
% Measure sampled extent and sheet area without inferring intrinsic PCA dimension.
centred = points-mean(points,2);
gram = centred.'*centred;
squareNorms = diag(gram);
distanceSquared = squareNorms+squareNorms.'-2*gram;
diameter = sqrt(max(0,max(distanceSquared(:))));
widths = (svd(centred.'/sqrt(size(points,2)))).';
edge1 = points(:,triangles(:,2))-points(:,triangles(:,1));
edge2 = points(:,triangles(:,3))-points(:,triangles(:,1));
normals = cross(edge1,edge2,1);
area = sum(sqrt(sum(normals.^2,1)))/2;
end

function [labels,representatives] = endpoint_clusters(joint,tolerance)
% Count tolerance clusters of PAIRS, not individual alternating endpoints.
labels = zeros(1,size(joint,2));
representatives = zeros(size(joint,1),size(joint,2));
numberOfClusters = 0;
for k = 1:size(joint,2)
    distances = sqrt(sum((representatives(:,1:numberOfClusters)-joint(:,k)).^2,1));
    match = find(distances <= tolerance,1);
    if isempty(match)
        numberOfClusters = numberOfClusters+1;
        representatives(:,numberOfClusters) = joint(:,k);
        match = numberOfClusters;
    end
    labels(k) = match;
end
representatives = representatives(:,1:numberOfClusters);
end

function family = directional_gap_family(bodies,gapGeometry,numberOfAngles)
% Vary a normal-section direction as an extra MODEL parameter at the fixed gap.
n = gapGeometry.normal;
T = [1,0;0,1;0,0];
assert(norm(T.'*n) < 1e-14,'This geometry uses the world-z gap normal.');
theta = linspace(0,pi,numberOfAngles);
directions = [cos(theta);sin(theta)];
radii = zeros(2,numberOfAngles);
for k = 1:2
    body = bodies(k);
    A = body.Q*diag(1./body.axes.^2)*body.Q.';
    gradient = A*(gapGeometry.contacts(:,k)-body.c);
    K = T.'*A*T/norm(gradient);
    normalCurvature = sum(directions.*(K*directions),1);
    radii(k,:) = 1./normalCurvature;
end
pair = sphere_gap_pair(radii(1,:),radii(2,:),gapGeometry);
family = struct('theta',theta,'radii',radii,'body1',pair.body1, ...
    'body2',pair.body2,'interpretation', ...
    'Directional models at the gap, NOT limits from different seeds of one map.');
end

function report_results(results)
% Distinguish established fixed-pair geometry from finite-resolution evidence.
fprintf('\n3D moving tangent-sphere Kelvin experiment\n');
fprintf('Filled focal disk: %d samples (interior rings + centre + boundary).\n', ...
    size(results.seeds,2));
fprintf('Exact supporting-plane gap %.6g; orbit tolerance %.1e.\n', ...
    results.gap.distance,results.options.orbitTol);
fprintf('One cycle = body 2 -> body 1 -> body 2.\n');
for k = 1:numel(results.runs)
    run = results.runs{k};
    fprintf('\nRule %s: %d cycles, converged=%d, pair clusters=%d (tol %.1e).\n', ...
        run.rule,run.cycles(end),run.converged, ...
        size(run.cluster_representatives,2),results.options.endpointTol);
    fprintf('  Final diameters [body1 body2]: %.3e %.3e\n',run.diameters(end,:));
    fprintf('  %s\n',run.fixed_set_status);
    fprintf('  Maximum true pair residual %.3e; reference errors %.3e %.3e\n', ...
        max(run.pair_residuals),run.reference_errors);
    fprintf('  Gap reference p1 = [%.10g %.10g %.10g]\n',run.reference.body1);
    fprintf('                p2 = [%.10g %.10g %.10g]\n',run.reference.body2);
    fprintf('  Full moving-map spectral radius %.8f; FD refinement difference %.3e\n', ...
        run.reference.linearisation.spectral_radius, ...
        run.reference.linearisation.step_refinement_difference);
    fprintf('  Image-domain violations (including repeated seeds): %d %d\n', ...
        sum(run.outside_counts,1));
    for j = 1:numel(results.options.collapseTols)
        fprintf('  First D <= %.1e: cycles [%.0f %.0f] (NaN = not reached).\n', ...
            results.options.collapseTols(j),run.collapse_cycles(j,:));
    end
    if any(run.outside_counts(:)) || ~all(run.reference.inside)
        warning('moving_kelvin:ImageOutside', ...
            'Rule %s: tangent-sphere inversion is not a domain-preserving reflection.',run.rule);
    end
    if ~run.converged
        warning('moving_kelvin:Unresolved', ...
            'Rule %s has not established a converged cloud at the requested tolerance.',run.rule);
    end
end
fprintf('\nPCA widths describe GLOBAL shape; a bent sheet can have 3 nonzero widths.\n');
fprintf('A small residual and diameter establish only tolerance-level point collapse.\n');
fprintf('An interior fixed pair is unique for each fixed local radius rule.\n');
fprintf('Directional-family intervals reflect model ambiguity, not seed dependence.\n');
end

function plot_results(results)
% Show material image sheets/curves, tolerance collapse, and sphere-rule sensitivity.
run = results.runs{1};
mesh = results.mesh;
colours = lines(numel(results.runs));
figure('Name','3D advected focal disk and image curves','Color','w');
layout = tiledlayout(2,2,'TileSpacing','compact');
for parity = 1:2
    ax = nexttile(layout,parity);
    hold(ax,'on');
    plot_bodies(ax,results.bodies);
    history = run.body1_images;
    if parity == 2
        history = run.body2_images;
    end
    shownCycles = unique([results.options.cyclesToShow,run.cycles(end)]);
    shownCycles = shownCycles(shownCycles <= run.cycles(end));
    if parity == 1
        shownCycles = shownCycles(shownCycles > 0);
    end
    cycleColours = turbo(max(2,numel(shownCycles)));
    for j = 1:numel(shownCycles)
        cycle = shownCycles(j);
        points = history(:,:,cycle+1);
        trisurf(mesh.triangles,points(1,:),points(2,:),points(3,:), ...
            'Parent',ax,'FaceColor',cycleColours(j,:), ...
            'FaceAlpha',0.18,'EdgeColor','none','HandleVisibility','off');
        % Every ring is mapped pointwise with the same mesh node indices.
        ringsToShow = unique([ceil(numel(mesh.ring_indices)/2),numel(mesh.ring_indices)]);
        for ring = ringsToShow
            indices = mesh.ring_indices{ring};
            closedIndices = indices([1:end,1]);
            visibility = 'off';
            if ring == ringsToShow(end)
                visibility = 'on';
            end
            plot3(ax,points(1,closedIndices),points(2,closedIndices),points(3,closedIndices), ...
                '-','Color',cycleColours(j,:),'LineWidth',1.2, ...
                'HandleVisibility',visibility,'DisplayName',sprintf('cycle %d',cycle));
        end
    end
    endpoint = run.reference.body1;
    if parity == 2
        endpoint = run.reference.body2;
    end
    plot3(ax,endpoint(1),endpoint(2),endpoint(3),'kp','MarkerFaceColor','k', ...
        'MarkerSize',10,'DisplayName','gap fixed-pair prediction');
    axis(ax,'equal'); grid(ax,'on'); view(ax,35,22);
    xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
    title(ax,sprintf('Mean-curvature images in body %d',parity));
    legend(ax,'Location','best');
end
ax = nexttile(layout,3);
hold(ax,'on');
for k = 1:numel(results.runs)
    current = results.runs{k};
    semilogy(ax,current.cycles,max(current.diameters(:,2),realmin), ...
        'Color',colours(k,:),'LineWidth',1.5,'DisplayName',current.rule);
end
for threshold = results.options.collapseTols(:).'
    yline(ax,threshold,':','HandleVisibility','off');
end
set(ax,'YScale','log'); grid(ax,'on');
xlabel(ax,'full cycles'); ylabel(ax,'sampled body-2 cloud diameter');
title(ax,'Point collapse is tolerance-dependent'); legend(ax,'Location','best');
ax = nexttile(layout,4);
hold(ax,'on');
styles = {'-','--',':'};
for j = 1:3
    semilogy(ax,run.cycles,max(run.rms_widths(:,j,2),realmin), ...
        styles{j},'LineWidth',1.4,'DisplayName',sprintf('RMS width %d',j));
end
set(ax,'YScale','log'); grid(ax,'on');
xlabel(ax,'full cycles'); ylabel(ax,'principal RMS width');
title(ax,'Global PCA is not intrinsic dimension'); legend(ax,'Location','best');
title(layout,'Advected focal disk: moving contacts and spheres, mean rule');

figure('Name','3D fixed-pair clouds and curvature ambiguity','Color','w');
layout = tiledlayout(2,2,'TileSpacing','compact');
ax = nexttile(layout,1);
hold(ax,'on');
family = results.directional_family;
plot3(ax,family.body1(1,:),family.body1(2,:),family.body1(3,:), ...
    'k-','LineWidth',3,'DisplayName','directional family (different models)');
plot3(ax,family.body2(1,:),family.body2(2,:),family.body2(3,:), ...
    'k-','LineWidth',3,'HandleVisibility','off');
for k = 1:numel(results.runs)
    current = results.runs{k};
    p1 = current.body1_images(:,:,end);
    p2 = current.body2_images(:,:,end);
    scatter3(ax,p1(1,:),p1(2,:),p1(3,:),20,colours(k,:),'filled', ...
        'DisplayName',sprintf('%s endpoints',current.rule));
    scatter3(ax,p2(1,:),p2(2,:),p2(3,:),20,colours(k,:),'filled', ...
        'HandleVisibility','off');
end
axis(ax,'equal'); grid(ax,'on'); view(ax,35,22);
xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
title(ax,'Do not pool models to infer one limit set'); legend(ax,'Location','best');
ax = nexttile(layout,2);
hold(ax,'on');
for k = 1:numel(results.runs)
    current = results.runs{k};
    semilogy(ax,1:current.cycles(end),max(current.cycle_changes,realmin), ...
        'Color',colours(k,:),'LineWidth',1.5,'DisplayName',current.rule);
end
yline(ax,results.options.orbitTol,'k:','DisplayName','orbit tolerance');
set(ax,'YScale','log'); grid(ax,'on');
xlabel(ax,'full cycles'); ylabel(ax,'largest parity cycle change');
title(ax,'Convergence of the actual moving map'); legend(ax,'Location','best');
ax = nexttile(layout,3);
hold(ax,'on');
for parity = 1:2
    for j = 1:2
        plot(ax,1:run.cycles(end),squeeze(run.representative_curvatures(j,parity,:)), ...
            'LineWidth',1.2,'DisplayName',sprintf('body %d, k_%d',parity,j));
    end
end
grid(ax,'on'); xlabel(ax,'full cycles'); ylabel(ax,'positive principal curvature');
title(ax,'Two curvatures at a moving contact (one seed)'); legend(ax,'Location','best');
ax = nexttile(layout,4);
hold(ax,'on');
for parity = 1:2
    semilogy(ax,run.cycles,max(run.triangular_areas(:,parity),realmin), ...
        'LineWidth',1.5,'DisplayName',sprintf('body %d',parity));
end
set(ax,'YScale','log'); grid(ax,'on');
xlabel(ax,'full cycles'); ylabel(ax,'advected triangular-mesh area');
title(ax,'Finite-time sheets can shrink to point limits'); legend(ax,'Location','best');
title(layout,'Fixed curvature rule versus a family of directional models');
end

function plot_bodies(ax,bodies)
% Render the physical ellipsoids as faint context for the material images.
[x,y,z] = sphere(28);
colours = [0.2,0.4,0.8;0.8,0.4,0.2];
for k = 1:2
    points = bodies(k).c+bodies(k).Q* ...
        (bodies(k).axes.*[x(:).';y(:).';z(:).']);
    surf(ax,reshape(points(1,:),size(x)),reshape(points(2,:),size(x)), ...
        reshape(points(3,:),size(x)),'FaceColor',colours(k,:), ...
        'FaceAlpha',0.07,'EdgeColor','none','HandleVisibility','off');
end
end

function checks = self_test
% Exercise projection, curvature, inversion, fixed-pair algebra, and metrics only.
[bodies,gapGeometry] = make_geometry(0.08);
options = struct('nTheta',8,'nRadial',2);
[seeds,mesh] = make_focal_seeds(bodies(2),options);
assert(isequal(size(seeds),[3,25]));
assert(all(mesh.radial_coordinate <= 1+10*eps));
query = seeds(:,[1,7,25]);
[contact,normal,curvatures,projectionError] = closest_jet(query,bodies(1));
curvatureError = 0;
for k = 1:size(query,2)
    direct = curvature_at_contact(contact(:,k),bodies(1));
    curvatureError = max(curvatureError,norm(curvatures(:,k)-direct));
end
assert(curvatureError < 2e-12);
assert(all(sum((query-contact).*normal,1) > 0));
rules = {'mean','min','max'};
gapResiduals = zeros(1,numel(rules));
for k = 1:numel(rules)
    [image,info] = moving_step(query,bodies(1),rules{k});
    offset = query-info.centre;
    directImage = info.centre+offset.* ...
        (info.radius.^2./sum(offset.^2,1));
    assert(max(abs(image(:)-directImage(:))) < 3e-12);
    reference = gap_reference(bodies,gapGeometry,rules{k});
    gapResiduals(k) = reference.residual;
    assert(reference.residual < 2e-12 && all(reference.inside));
    assert(reference.linearisation.step_refinement_difference < 2e-7);
    assert(reference.linearisation.spectral_radius < 1);
end

% A three-node, three-cycle component test exercises history storage and
% collapse diagnostics, not the default disk experiment or any plotting.
tinyOptions = struct('maxCycles',3,'orbitTol',1e-12,'endpointTol',1e-7, ...
    'collapseTols',[1e-3,1e-6]);
tinyMesh = struct('triangles',[1,2,3]);
reference = gap_reference(bodies,gapGeometry,'mean');
tinyRun = iterate_cloud(query,tinyMesh,bodies,'mean',reference,tinyOptions);
assert(isequal(size(tinyRun.body1_images),[3,3,4]));
assert(isequal(size(tinyRun.body2_images),[3,3,4]));
assert(isequal(size(tinyRun.rms_widths),[4,3,2]));
assert(all(isfinite(tinyRun.diameters(2:end,:)),'all'));
assert(all(isfinite(tinyRun.pair_residuals)));
assert(all(tinyRun.diameters(end,:) < tinyRun.diameters(2,:)));
assert(~tinyRun.converged);

% A sphere has equal curvatures everywhere and the tangent sphere is itself.
sphereBody = struct('axes',2*ones(3,1),'Q',eye(3),'c',[0.2;-0.1;0.4]);
sphereQueries = sphereBody.c+[3,0,1;0,-4,2;0,0,3];
[sphereImages,sphereInfo] = moving_step(sphereQueries,sphereBody,'mean');
offset = sphereQueries-sphereBody.c;
expected = sphereBody.c+4*offset./sum(offset.^2,1);
assert(max(abs(sphereImages(:)-expected(:))) < 2e-12);
assert(max(abs(sphereInfo.curvatures(:)-0.5)) < 2e-8);

% Distinguish planar sample width from a sheet bent in three world directions.
[diameter,widths,area] = cloud_metrics([0,1,0;0,0,1;0,0,0],[1,2,3]);
assert(abs(diameter-sqrt(2)) < 20*eps && widths(3) < 20*eps);
assert(abs(area-0.5) < 20*eps);
paraboloidNodes = [0,1,-1,0,0;0,0,0,1,-1;0,1,1,1,1];
[~,bentWidths] = cloud_metrics(paraboloidNodes, ...
    [1,2,4;1,4,3;1,3,5;1,5,2]);
assert(all(bentWidths > 0.1),'A bent 2D sheet need not have zero third PCA width.');
family = directional_gap_family(bodies,gapGeometry,9);
assert(max(abs(family.body1(1:2,:)),[],'all') < 20*eps);
assert(max(abs(family.body2(1:2,:)),[],'all') < 20*eps);
checks = struct('passed',true,'sample_count',size(seeds,2), ...
    'projection_errors',projectionError,'curvature_error',curvatureError, ...
    'gap_pair_residuals',gapResiduals,'experiment_executed',false);
fprintf('Moving-sphere 3D self_test passed: 25 focal nodes, 3 test queries, 3 rules.\n');
fprintf('Three-node, three-cycle storage/diagnostic component test passed.\n');
fprintf('Projection errors %.3e %.3e; curvature error %.3e; gap residual %.3e.\n', ...
    projectionError,curvatureError,max(gapResiduals));
end
