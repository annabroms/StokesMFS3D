function [mode_map,info] = papavassiliou_squirmer_passive_map(h,a,series_opts)
%PAPAVASSILIOU_SQUIRMER_PASSIVE_MAP Exact coaxial unit-sphere velocity map.
%
%   [MODE_MAP,INFO] = PAPAVASSILIOU_SQUIRMER_PASSIVE_MAP(H,A,SERIES_OPTS)
%   returns the exact bispherical-series response of a spherical squirmer
%   interacting coaxially with a passive sphere. This implementation is
%   restricted to two equal unit spheres: A must be a finite scalar equal
%   to one within roundoff, and H must be a finite positive scalar giving
%   the surface-to-surface gap.
%
%   The particle centres and positive-axis convention are
%
%       c^(1) = -(1 + H/2) e_z,   c^(2) = (1 + H/2) e_z.
%
%   Particle 1 is the squirmer, particle 2 is passive, and positive
%   velocity points from particle 1 toward particle 2. MODE_MAP is 2-by-2
%   and is ordered so that
%
%       [V_sq; V_ps] = MODE_MAP * [B1; B2].
%
%   Here B1 and B2 are the tangential squirming coefficients in equations
%   (4.2) and (B3)--(B4) of Papavassiliou and Alexander,
%   "Exact solutions for hydrodynamic interactions of two squirming
%   spheres", J. Fluid Mech. 813 (2017), 618--646. For alpha=0 their
%   physical surface-distortion velocity on particle 1 is
%
%       u_s = [B1 + B2 (e_z dot n)] sin(theta) e_theta.
%
%   The physical slip array passed to solve_mobility_with_DLP is therefore
%
%       u_slip = u_s
%              = -[B1 + B2 (e_z dot n)] (I - n*n.') e_z,
%
%   and the repository boundary condition is
%
%       u = V + Omega x (x-c) + u_slip.
%
%   The calculation uses the translating-sphere coefficient equations
%   (3.12)--(3.13), the conjugate forces (3.16), the reciprocal theorem
%   (4.1), and the tangential-mode integrals (B3)--(B4). Since alpha=0,
%   P1(cos(alpha)) = P2(cos(alpha)) = 1.
%
%   SERIES_OPTS is an optional scalar struct with fields:
%       tol               relative term tolerance in (0,1) (default 1e-13)
%       min_terms         positive integer minimum terms (default 20)
%       consecutive_terms positive integer consecutive small terms
%                         required for convergence (default 8)
%       max_terms         positive integer maximum terms (default 2000)
%
%   The resistance and reciprocal-integral series are accumulated
%   together. Convergence is declared only after the largest normalized
%   contribution has remained below TOL for CONSECUTIVE_TERMS terms after
%   MIN_TERMS. The function errors with identifier
%   papavassiliou_squirmer_passive_map:seriesNotConverged if MAX_TERMS is
%   reached first.
%
%   INFO contains:
%       geometry             H, A, xi1, xi2, R, and the two centres
%       series_opts          validated adaptive-series options
%       terms                number of retained terms
%       converged            true on successful return
%       tail_estimate        largest normalized final term
%       resistance_matrix    force map for the two conjugate translations
%       reciprocal_rhs       Appendix-B integrals for unit B1 and unit B2
%       reciprocity_defect   relative symmetry defect of the resistance
%       condition_number     2-norm condition number of the resistance
%
%   Example:
%       h = 0.2;
%       B1 = 3/2;
%       beta = 0;
%       [mode_map,info] = papavassiliou_squirmer_passive_map(h,1);
%       velocities = mode_map * [B1; beta*B1]
%       assert(info.converged)
%
%   In the far-field limit,
%
%       V_sq -> 2*B1/3,   V_ps -> 0.
%
%   Calling PAPAVASSILIOU_SQUIRMER_PASSIVE_MAP with no arguments runs a
%   deterministic self-test of convergence, resistance reciprocity,
%   B1/B2 linearity, tolerance stability, and the far-field limit.
%
%   See also SOLVE_MOBILITY_WITH_DLP.

if nargin == 0
    self_test();
    mode_map = [];
    info = [];
    return;
end
if nargin < 2 || isempty(a)
    a = 1;
end
if nargin < 3 || isempty(series_opts)
    series_opts = struct();
end

validateattributes(h,{'numeric'}, ...
    {'real','finite','scalar','positive'},mfilename,'h',1);
validateattributes(a,{'numeric'}, ...
    {'real','finite','scalar','positive'},mfilename,'a',2);
if abs(a-1) > 32*eps(max(1,abs(a)))
    error('papavassiliou_squirmer_passive_map:unitSphereOnly', ...
        ['Only two unit spheres are supported. The analytical-radius ', ...
         'input a must equal one; received a = %.17g.'],a);
end
if ~isstruct(series_opts) || ~isscalar(series_opts)
    error('papavassiliou_squirmer_passive_map:badSeriesOptions', ...
        'series_opts must be a scalar struct.');
end
series_opts = fill_series_defaults(series_opts);

xi1 = acosh(1+h/2);
xi2 = -xi1;
R = sinh(xi1);
centres = [0 0 -(1+h/2); 0 0 (1+h/2)];

resistance_sum = zeros(2);
reciprocal_sum = zeros(2);
small_term_count = 0;
tail_estimate = inf;
converged = false;
mu = 1;

force_prefactor = 4*pi*mu*sqrt(2*R);
B1_prefactor = 8*pi*mu*sqrt(2*R)/5;
B2_prefactor = 8*pi*mu*sqrt(2*R)/105;

for ell = 1:series_opts.max_terms
    [A,B,C,D] = translation_coefficients_equal_spheres(ell,xi1,R);
    ell_weight = ell*(ell+1);

    resistance_term = -force_prefactor*ell_weight * ...
        [A+B; C+D];

    exp_x = exp(xi1);
    exp_minus_x = exp(-xi1);
    exp_minus_2ell_x = exp(-2*ell*xi1);
    exp_minus_2ellp2_x = exp(-2*(ell+1)*xi1);

    AD = exp_x*A + exp_minus_2ellp2_x*D;
    BC = exp_minus_x*B + exp_minus_2ell_x*C;

    B1_bracket = ...
        ((4*ell+7)*sinh(xi1)-cosh(xi1))*AD + ...
        ((4*ell-3)*sinh(xi1)-cosh(xi1))*BC;

    B2_AD_factor = ...
        8*(2*ell^2+ell-3)*exp(2*xi1) - ...
        (32*ell^2+67*ell+27) + ...
        (ell+2)*(16*ell+27)*exp(-2*xi1);
    B2_BC_factor = ...
        (ell-1)*(16*ell-11)*exp(2*xi1) - ...
        (32*ell^2-3*ell-8) + ...
        8*(ell+2)*(2*ell-1)*exp(-2*xi1);
    B2_bracket = B2_AD_factor*AD + B2_BC_factor*BC;

    reciprocal_term = ell_weight * ...
        [B1_prefactor*B1_bracket(:), ...
         B2_prefactor*B2_bracket(:)];

    if ~all(isfinite(resistance_term(:))) || ...
            ~all(isfinite(reciprocal_term(:)))
        error('papavassiliou_squirmer_passive_map:nonfiniteSeries', ...
            ['The bispherical series produced a non-finite term at ', ...
             'ell = %d for h = %.17g.'],ell,h);
    end

    resistance_sum = resistance_sum + resistance_term;
    reciprocal_sum = reciprocal_sum + reciprocal_term;

    resistance_scale = max(abs(resistance_sum),[],2);
    reciprocal_scale = max(abs(reciprocal_sum),[],1);
    resistance_rel = max(abs(resistance_term) ./ ...
        max(resistance_scale,realmin),[],'all');
    reciprocal_rel = max(abs(reciprocal_term) ./ ...
        max(reciprocal_scale,realmin),[],'all');
    tail_estimate = max(resistance_rel,reciprocal_rel);

    if ell >= series_opts.min_terms && tail_estimate < series_opts.tol
        small_term_count = small_term_count + 1;
    else
        small_term_count = 0;
    end

    if small_term_count >= series_opts.consecutive_terms
        converged = true;
        break;
    end
end

if ~converged
    error('papavassiliou_squirmer_passive_map:seriesNotConverged', ...
        ['Bispherical series did not reach tol = %.3e after %d terms ', ...
         '(last normalized term %.3e, h = %.17g).'], ...
        series_opts.tol,series_opts.max_terms,tail_estimate,h);
end

% Convert the assembled Appendix-B modal columns to the paper's B1/B2 basis.
paper_mode_sign = diag([1 -1]);
mode_map = resistance_sum.' \ (reciprocal_sum * paper_mode_sign);




if ~all(isfinite(mode_map(:)))
    error('papavassiliou_squirmer_passive_map:singularResistance', ...
        'The transposed resistance solve produced a non-finite velocity map.');
end

info = struct();
info.geometry = struct( ...
    'h',h,'a',a,'xi1',xi1,'xi2',xi2,'R',R,'centres',centres);
info.series_opts = series_opts;
info.terms = ell;
info.converged = converged;
info.tail_estimate = tail_estimate;
info.resistance_matrix = resistance_sum;
info.reciprocal_rhs = reciprocal_sum;
info.reciprocity_defect = norm( ...
    resistance_sum-resistance_sum.','fro') / ...
    max(norm(resistance_sum,'fro'),realmin);
info.condition_number = cond(resistance_sum);

end

function opts = fill_series_defaults(opts)
defaults = struct( ...
    'tol',1e-13, ...
    'min_terms',20, ...
    'consecutive_terms',8, ...
    'max_terms',2000);
names = fieldnames(defaults);
for k = 1:numel(names)
    name = names{k};
    if ~isfield(opts,name) || isempty(opts.(name))
        opts.(name) = defaults.(name);
    end
end

allowed = fieldnames(defaults);
provided = fieldnames(opts);
unknown = setdiff(provided,allowed);
if ~isempty(unknown)
    error('papavassiliou_squirmer_passive_map:badSeriesOptions', ...
        'Unknown series_opts field "%s".',unknown{1});
end

validateattributes(opts.tol,{'numeric'}, ...
    {'real','finite','scalar','positive'},mfilename,'series_opts.tol');
if opts.tol >= 1
    error('papavassiliou_squirmer_passive_map:badSeriesOptions', ...
        'series_opts.tol must be smaller than one.');
end
validate_integer_option(opts.min_terms,'min_terms');
validate_integer_option(opts.consecutive_terms,'consecutive_terms');
validate_integer_option(opts.max_terms,'max_terms');
if opts.max_terms < opts.min_terms + opts.consecutive_terms - 1
    error('papavassiliou_squirmer_passive_map:badSeriesOptions', ...
        ['series_opts.max_terms must be at least min_terms + ', ...
         'consecutive_terms - 1.']);
end
end

function validate_integer_option(value,name)
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
        ~isfinite(value) || value < 1 || value ~= round(value)
    error('papavassiliou_squirmer_passive_map:badSeriesOptions', ...
        'series_opts.%s must be a positive integer scalar.',name);
end
end

function [A,B,C,D] = translation_coefficients_equal_spheres(ell,xi,R)
% Coefficients for conjugate translations (V1,V2)=(1,0),(0,1).
aexp = ell + 3/2;
bexp = ell - 1/2;
cexp = ell + 1/2;

plus_matrix = 2 * [ ...
    cosh(aexp*xi), cosh(bexp*xi); ...
    (2*ell+3)*sinh(aexp*xi), (2*ell-1)*sinh(bexp*xi)];
minus_matrix = 2 * [ ...
    sinh(aexp*xi), sinh(bexp*xi); ...
    (2*ell+3)*cosh(aexp*xi), (2*ell-1)*cosh(bexp*xi)];

f = exp(-aexp*xi)/(2*ell+3) - ...
    exp(-bexp*xi)/(2*ell-1);
rhs_base = (R/sqrt(2)) * [ ...
    f; ...
    2*exp(-cexp*xi)*sinh(xi)];

plus_coefficients = plus_matrix \ (rhs_base * [1 1]);
minus_coefficients = minus_matrix \ (rhs_base * [1 -1]);

A = 0.5*(plus_coefficients(1,:) + minus_coefficients(1,:));
D = 0.5*(plus_coefficients(1,:) - minus_coefficients(1,:));
B = 0.5*(plus_coefficients(2,:) + minus_coefficients(2,:));
C = 0.5*(plus_coefficients(2,:) - minus_coefficients(2,:));
end

function self_test()
opts = struct('tol',1e-12);
[map_near,info_near] = papavassiliou_squirmer_passive_map(0.2,1,opts);
assert(info_near.converged);
assert(info_near.reciprocity_defect < 1e-10);

coeff_a = [3/2; -0.4];
coeff_b = [-0.2; 0.7];
linearity_defect = norm( ...
    map_near*(coeff_a+coeff_b) - ...
    (map_near*coeff_a+map_near*coeff_b),inf);
assert(linearity_defect < 100*eps*max(1,norm(map_near,inf)));

[map_tight,info_tight] = papavassiliou_squirmer_passive_map( ...
    0.2,1,struct('tol',1e-13));
tolerance_defect = norm(map_near-map_tight,inf) / ...
    max(norm(map_tight,inf),realmin);
assert(info_tight.converged);
assert(tolerance_defect < 1e-9);

[map_far,info_far] = papavassiliou_squirmer_passive_map(1000,1);
far_response = map_far*[1;0];
assert(info_far.converged);
assert(abs(far_response(1)-2/3) < 5e-7);
assert(abs(far_response(2)) < 5e-7);

fprintf(['papavassiliou_squirmer_passive_map self-test passed: ', ...
    'terms(h=0.2)=%d, reciprocity=%.3e, stability=%.3e, ', ...
    'far response=[%.12g, %.3e].\n'], ...
    info_tight.terms,info_tight.reciprocity_defect,tolerance_defect, ...
    far_response(1),far_response(2));
end
