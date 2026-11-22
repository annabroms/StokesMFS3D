# Wall and close-pair image solvers

[solve_resistance_wall_images](../resistance/solve_resistance_wall_images.m)
computes forces/torques from prescribed rigid motion.
[solve_mobility_wall_images](../mobility/solve_mobility_wall_images.m)
computes rigid motion from prescribed forces/torques. Both support equal
spheres above a horizontal no-slip plane, with simultaneous close wall and
particle contacts, zero imposed slip and zero background flow.

Wall lines default to **ST_tD**. Close-pair lines default to **SDR**, with one interior
line owned by each member of a neighbouring pair. Every source, including
each pair source, includes its complete wall correction. Wall reflections
are prescribed by the Green function and are not additional fitted unknowns.

Change the two combinations independently in the same options structure
passed to either wall solver:

```matlab
opt.wallImageFamilies = "ST_tD";  % Default wall combination.
opt.pairImageFamilies = "SDR";    % Default particle-pair combination.
opt.viscosity = 1;               % Default dynamic viscosity.
opt.SVDsolve = false;            % Default: one-body Helsing GMRES.
```

For example, `"SRTD"` includes every family and all three stresslet
orientations; `"SD"` uses only line Stokeslets and dipoles; `"none"` disables
that type of line and its additional fitting caps. Standard proxy Stokeslets
remain present. Tokens can also be separated with `+`.

`T_t` is the axial stresslet direction along the relevant line. `T_n` and
`T_n2` are its two transverse directions, and `T` includes all three. On
a wall line t points toward the wall; on each particle-pair line t points
toward that neighbour. Source strengths remain three Cartesian components.
Each selected vector family adds three coefficients per line point. Broader
combinations may require a larger surface grid or capCount to retain the
documented oversampling. The plots follow the selected families automatically.
Combination strings select sets of families; coefficient ordering is
deterministic and is recorded by the basis labels. Default ordering is unchanged.

The separate free-space solvers are
[solve_resistance_images](../resistance/solve_resistance_images.m) and
[solve_mobility_images](../mobility/solve_mobility_images.m). Their existing
single selector is `opt.imageFamilies="SRD"` by default, an equivalent set
to SDR. Their `T` defaults to the axial stresslet (the older
use_orth_dir_images option changes it to one transverse direction); the wall solvers' explicit
T_t/T_n/T_n2 notation distinguishes the orientations.

## Global dense SVD option

Set `opt.SVDsolve=true` in either wall solver to assemble and solve the
global rectangular system using the relative truncated-SVD threshold opt.tol.
The default false uses one-body preconditioning and helsing_gmres. The dense
resistance system is the complete physical velocity matrix; mobility uses
the complete recompleted matrix with the same force/torque projection and
completion source as the iterative formulation. It has no identity-diagonal
replacement and no local preconditioner factorizations.

The dense path uses direct kernels even when opt.fmm=true, and therefore
needs no compiled FMM dependency. It prints matrix dimensions and a memory
estimate before assembly. Global storage and factorization can be expensive.
It returns iters=0, details.real_res=NaN and details.converged=NaN because
the GMRES stopping criterion does not apply. The relative rectangular
residual is details.linearSystemResidual; independently measured surface
errors are returned as usual, without accuracy assertions. Global SVD rank
is details.retainedRanks. Global assembly/factorization timings are reported
separately; solveSeconds includes both in dense mode.

`opt.viscosity` defaults to 1 and sets physical force/source units in both
modes. The former opt.mu is accepted as an alias when viscosity is omitted;
an explicit viscosity takes precedence. The free-space image solvers also
default to opt.viscosity=1 and scale resistance loads and mobility consistently.

```matlab
setpath;                         % StokesHalfspaceFMM must be on the path.
blake.setup;                     % Select the installed compiled FMM3D.
q = [0,0,1.001; 2.001,0,1.001; 1,4,1.002];
[rin,rout] = init_spheres(q,0.63,700,2);
opt = struct('fmm',1,'N_image',30,'tol',1e-12,'gmres_tol',1e-9);
U = repmat([1;0;-1;0;0.2;0],3,1);
[Fvec,iters,lambda_norm,err_res,lambda,details] = ...
    solve_resistance_wall_images(q,rin,rout,U,opt);

Fvec = repmat([1;0;-1;0;0;0.2],3,1);
[U,iters,lambda_norm,uerr,lambda,details] = ...
    solve_mobility_wall_images(q,rin,rout,Fvec,opt);
```

The general [init_spheres](../geometry/spheres/init_spheres.m) replaces the
former wall-specific initializer. Its optional fifth input is physical
radius: `init_spheres(q,Rp,N,a_glob,a)`. Rp is a fraction of a; default a=1.
Wall checks belong to solver setup. Old four-input calls retain their
original unit-radius behavior. If using the initializer's third output,
opt.radius records a and its legacy default is `opt.fmm=0`. Set `opt.fmm=1`
to select the hybrid backend, and choose the solver tolerances explicitly.

Call each solver without arguments to run its own deterministic self-test:

```matlab
solve_resistance_wall_images;
solve_mobility_wall_images;
```

No-argument calls also return their self-test results when outputs are
requested, making the reported grids, coefficients and residuals inspectable.

Both tests place three spheres close to the wall; the first two also have
particle surface gap 1e-3. They leave three figures open: discretization
with wall and pair close-ups, all seven source-type magnitudes, and surface
residual at every collocation node. The mobility test uses fixed forces and
torques independently of the resistance result. Missing the 1e-6 surface
accuracy goal is printed as an accuracy result, without failing the test.
Ordinary calls plot only with `opt.plot=1`.

## Coordinates, coefficients and diagnostics

`rin` and `rout` supplied to the solvers contain ordinary proxy and surface
grids, grouped by particle. Wall/pair lines and contact caps are appended
internally. The final physical coordinates are `details.rin` and
`details.rout`. Shared family positions appear more than once in rin because
each family has its own coefficients.

N is the initial proxy point count per particle. `N_image` counts geometric
points per line. Each default line adds three vector families, hence
9*N_image scalar unknowns. In general a line adds
3*N_image*(number of selected vector families) scalar unknowns; T selects
three such families. `N_total(k)` counts all source vectors on particle k;
`M(k)` counts its final surface collocation nodes. These counts are returned
in details and can differ across particles. The fit must have M(k)>=2*N_total(k).
Increase the supplied surface grid or capCount if setup reports insufficient
oversampling; no grid or SVD threshold is silently changed.

`U` and `Fvec` have six entries per particle, respectively
`[vx;vy;vz;omega_x;omega_y;omega_z]` and `[Fx;Fy;Fz;Tx;Ty;Tz]`, with
rotations/torques about the particle centre. Fvec is force/torque **on the
fluid**, so drag on the particle is -Fvec.

The returned lambda is the **total physical** mixed coefficient vector.
Components are interleaved xyz by node and grouped by body, then by family.
The default order is:

```text
proxy S; wall S; wall T_t; wall D;
pair 1 S; pair 1 D; pair 1 R; pair 2 S; pair 2 D; pair 2 R; ...
```

Inactive blocks are omitted and neighbours are sorted by particle index.
`details.problem.body{k}.basis` records labels, locations and source types
in the actual selected order, with normalized coordinates. On wall lines
T_t fixes the direction toward the wall, (0,0,-1), while fitting all three
Cartesian strengths. On pair lines it follows that neighbour's direction.
Package calls receive -lambda for T_t and -2*lambda for D to match the
existing MFS convention.

Internally coordinates relative to the wall are divided by radius a,
viscosity is one, and rigid data are `[V;a*Omega]`. Physical S strengths
equal mu*a times the dimensionless strengths, T_t/R equal mu*a^2 times,
and D equal mu*a^3 times. `details.sourceUnits` records this conversion per
scalar coefficient; `details.lambda` stores the dimensionless coefficients.
`lambda_norm` is their scalar infinity norm. The source figures instead
show vector magnitudes, separately for proxy S, wall S/T_t/D and pair S/D/R.

For example, evaluate the final collocation flow using prepared geometry:

```matlab
uvec = apply_wall_images_velocity( ...
    lambda./details.sourceUnits,details.problem);
```

Velocity output is interleaved by collocation node. A custom evaluation
grid passed to this shared evaluator must be in normalized coordinates.

`err_res`/`uerr` is the maximum independent pointwise velocity discrepancy
divided by the maximum sampled rigid speed. The independent grid includes
full-sphere quadrature, broad and inner Fibonacci caps toward every wall
and close neighbour, and closest poles. Zero rigid speed selects absolute
velocity error. `compute_residual=0` returns NaN for this independent error;
the collocation check still runs. Diagnostics use direct kernels independently
of FMM. `details.collocationResidual` aligns exactly with final rout.
`details.real_res` is the separately recomputed square GMRES residual.
A small GMRES residual does not certify resolved surface velocity or forces.

## Preconditioner and matvec

In the default one-body mode, each body is factored independently, including all of its wall/pair sources
and its wall and pair collocation caps. No two-body preconditioner is used.
The SVD factors retain modes strictly above `opt.tol` times the largest
singular value. Y and UU use the same convention as the main solvers:
`Y{k}*(UU{k}*x)` maps a body's collocation unknown to its source coefficients.
UU stores retained left singular vectors in rows; Y stores retained right
vectors divided by singular values. Kin' maps mixed coefficients to force
and torque; Kout maps rigid motion to surface velocity. Rotlets contribute
torque, and T_t/D contribute no force or torque. Mobility applies the
force/torque projection using six-column maps, without a full LL matrix.
Both solvers use helsing_gmres. They enable its optional second
orthogonalization pass and check the recomputed residual before early
termination, because close-pair systems can lose Krylov orthogonality.
If the recurrence has reached tolerance but the true residual stops improving
for five checks, helsing_gmres returns the best checked iterate and reports
its true residual; the solver's converged flag remains false when appropriate.

The production matvec uses one global half-space S action and one blocked
action for each active optional family. When fmm=1, S is accelerated;
T_t, D and R use the installed package's **direct-only** kernels. This is
a hybrid backend, with quadratic optional-family work at fixed per-body
resolution and neighbour count. No global dense interaction matrix or
particle-pair evaluation loop is built into the matvec. With fmm=0 all
families use direct kernels. Requesting FMM with an unavailable compiled
dependency raises the dependency's error; there is no silent fallback.

The self correction uses retained SVD vectors, rather than dense local
velocity matrices during iteration. Identity replacement of rectangular
blocks can have regularization bias; independently inspect surface errors.
Prepared geometry does not imply reusable FMM trees. The diagnostic problem
retains the local factors and maps and can therefore consume substantial
memory at large per-body resolutions.

At very small gaps, retained near-null modes can amplify source cancellation
and evaluation roundoff. Tightening the FMM tolerance alone may then leave
the action error unchanged. Compare the true GMRES residual with its recurrence
estimate and inspect the source magnitudes; a truncated SVD cutoff is part of
the numerical method, not just an implementation setting. Both `converged`
and `real_res` are reported, without asserting a surface accuracy target.

The full flat option list, defaults, dimensions and normalization are in
[prepare_wall_images_problem](prepare_wall_images_problem.m). Timings for
assembly, factorization, preparation, solve, validation and a warm matvec
are returned in `details.timing`. `gmres_tol`, `maxit`, `fmm_tol`,
`gmres_verbose`, `plot` and `compute_residual` follow main-solver naming.

Run the small focused algebraic checks with:

```matlab
report = test_wall_images_solvers;
report = test_wall_images_svd;         % Global dense path and physical scaling.
report = test_image_solver_viscosity;  % Renamed free-space solver APIs.
```

They compare against independent dense square systems, test mixed-family
moments, direct/hybrid agreement, zero loads and physical scaling. Their
deliberately coarse grids do not assert a 1e-6 surface error.

A full-size action check on the self-test geometry (2,508 coefficient vectors,
8,564 targets, unit viscosity, FMM tolerance 1e-12) used the deterministic
coefficient vector `cos((1:numel(details.lambda)).')`. Relative hybrid/direct
velocity discrepancy was 3.965e-14. Warm calls measured 1.32 s hybrid and
1.72 s direct in one run; these timings include all source families but not
local SVD factor application. They are individual measurements, not a scaling
benchmark. The self-test separately reports the complete preconditioned
matvec time and actual solve accuracy.

An earlier three-sphere validation with N_image=30, tol=1e-12 and both gaps
equal to 1e-3 gave the following one-body iterative results:

| Self-test | Iterations | True GMRES residual | Independent surface error |
| --- | ---: | ---: | ---: |
| Resistance | 260 | 2.200e-3 | 1.156e-2 |
| Mobility | 249 | 1.543e-6 | 9.399e-3 |

Both converged flags were false at requested GMRES tolerance 1e-9. The
mobility force/torque constraint error was 5.457e-12. These are deliberately
difficult simultaneous gaps of 1e-3, and the current defaults are **unresolved**
against the 1e-6 surface goal. A successful algebraic regression does not
certify this physical accuracy. The self-tests keep their figures and print
the measured errors without treating the surface goal as a test assertion.

Saved interactive figure bundles and PNG exports from these runs are in
[visuals/](../visuals/) with the `wall_images_resistance_` and
`wall_images_mobility_` prefixes. Suffixes 1, 2 and 3 denote discretization,
source magnitudes and collocation residual, respectively.
