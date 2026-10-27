# Proxy-only Stokes flow above a planar wall

The wall solvers impose rigid motion on equal spheres above an infinite,
stationary, no-slip plane. They use standard interior proxy Stokeslets and
the complete Blake Green function. No contact image lines or enhanced
source families are added. Background flow and prescribed slip are zero.

The half-space evaluation is provided by the independent sibling repository
`../StokesHalfspaceFMM`. Its fast representation is due to **Z. Gimbutas,
L. Greengard and S. Veerapaneni**, *Simple and efficient representations for
the fundamental solutions of Stokes flow in a half-space*, JFM 776 (2015), R1,
<https://doi.org/10.1017/jfm.2015.302>. Its independent direct reference uses
J. R. Blake (1971), equation (16), <https://doi.org/10.1017/S0305004100049902>.

## Demos

Run the demos from the StokesMFS3D root. Both use `setpath`, which adds the
sibling `../StokesHalfspaceFMM` MATLAB package.

`demo/demo_wall_mobility.m` is the small direct-Blake example. It solves
mobility for five equal-height spheres with random force/torque data, then
solves resistance using the resulting velocities. It reports both independent
surface residuals and the dimensionless force round-trip error. Run it with:

```matlab
run('demo/demo_wall_mobility.m')
```

`demo/demo_random_wall_torque.m` defaults to 100 randomly packed spheres in a
layer with center-height span 5, unit z-torque on every sphere, minimum body
gap 1, and wall gap 0.5. The full vertical extent of the spheres is therefore
7. It visualizes the geometry before solving, including the wall. After the
mobility solve, a three-panel figure colors the sphere surfaces by rigid-body
speed, scatters the independent relative surface residual, and shows the
magnitude of Blake traction using a separate `hot` colormap; all panels show
the wall. Pass a numeric second argument
to set the center-height span; the 2D `grow_cluster` layout preserves the
minimum body gap. For example, `demo_random_wall_torque(100,2)` gives center
heights spanning 2 and a full sphere extent of 4. Pass `[]` explicitly for a
3D cloud. The default FMM backend requires compiled FMM3D; initialize it before
running:

```matlab
setpath;
blake.setup;
results = demo_random_wall_torque;
```

For a quick small run, use `demo_random_wall_torque(5)`.

The fourth solver output is the independent surface residual. A small
round-trip error does not guarantee a small surface residual.

The Blake action backend is selected with `backend`: `"blake"` evaluates
Blake's formula directly, `"ggv-direct"` uses the independent direct half-space
representation, and `"fmm"` uses the accelerated half-space evaluator. The
`"fmm"` backend requires compiled FMM3D, initialized with `blake.setup`, and
never falls back to a direct backend. The preconditioner is selected separately
with `preconditioner="wall"` or `"freeSpace"`. Both use the full Blake
interaction in the physical solve; `"freeSpace"` omits the wall only from the
one-body preconditioner.

The first five outputs follow the basic solvers: answer, iteration count,
physical source infinity norm, independent boundary error, and total physical
sources. The sixth output contains convergence, spectra/ranks, timing, actual
gaps and validation samples. The basic solvers also accept these options with
`kernel="blake"`; ordinary free-space calls retain their original behavior.
The wall mobility wrapper has no slip argument. Explicit Blake selection in
the basic mobility solver requires `u_slip=[]`.

Coordinates are particle-grouped `NP-by-3` or `MP-by-3`. Vector components are
interleaved by node: `[x1;y1;z1;x2;y2;z2;...]`. Rigid motion uses six entries
`[Vx;Vy;Vz;Omega_x;Omega_y;Omega_z]` per particle, and forces use
`[Fx;Fy;Fz;Tx;Ty;Tz]`. Rotations and torques are about the particle centre.

`FT` is the force/torque **on the fluid**. At force balance it equals the
external load applied to the particle; hydrodynamic drag on the particle is
`-FT`. Each source strength is a point force, without a quadrature weight:

```text
F = sum_n lambda_n
T = sum_n (source_n-centre) cross lambda_n
```

All public arrays are dimensional. Internally, coordinates relative to the wall
are divided by radius `a` and viscosity is one. Taking velocity unit one,
internal velocities are `[V;a*Omega]`, internal forces are
`[F/(mu*a);T/(mu*a^2)]`, and physical strengths equal `mu*a` times internal
strengths. The velocity kernel includes `1/(8*pi*mu)` exactly once.

Only equal spheres are supported: radius may be a scalar or an equal-valued
vector of length P. Particles must have positive body gaps and wall clearance.

## One-body preconditioning and dense references

Choose `preconditioner="wall"` (default) or `"freeSpace"`. Both retain the
complete Blake kernel in the physical global interaction. The free-space
option omits the wall only from the isolated-body preconditioner.

With `H=K_in'`, `C=H'/(H*H')` and `Q=I-C*H`, the resistance block is
`B=Gpre`, and the mobility block is `B=Gpre*Q+K_out*H`, both 3M-by-3N.
The pseudoinverse keeps singular values strictly greater than
`svdCutoff*max(singularValues)`. Defaults are **svdCutoff=1e-12** and
**proxyRadius=0.63**; proxy radius is a fraction of the physical sphere radius.
Retained spectra and ranks are returned. No extra rank-one correction,
symmetrization, diagonal shift, or singular-kernel regularization is added.

The dedicated resistance and mobility wall solvers use `helsing_gmres` on
target-space unknowns x. Its callback is the shared matrix-free operator; each
iteration applies the global physical interaction through the selected Blake
backend and corrects the per-body self block without assembling a global dense
matrix. For auxiliary sources `z=B^+*x`, the square operator is
`x+Aphysical*z-diag(B)*z`. Mobility's `K_out*H` terms cancel in the callback.
The free-space preconditioner still retains the physical wall self correction.

Resistance reconstructs `lambda=z`, `FT=H*lambda`. Mobility constructs
`lambda0=C*FT`, solves against `-GBlake*lambda0`, then reconstructs
`lambda=Q*z+lambda0`, `U=-H*z`. Total physical strengths are always returned.

Identity replacement regularizes a rectangular self problem. Thus the two
preconditioner choices can have different discretization/regularization bias.
Compare their physical errors before interpreting iteration counts.
`solverMode="denseSVD"` computes a global rectangular truncated SVD without
identity replacement, and is restricted to P<=2. Dense mobility uses the
recompleted rectangular operator, not an identity-diagonal approximation.

For repeated loads on the same geometry:

```matlab
opt.preparedProblem = prepare_wall_problem(q,rin,rout,opt);
```

Pass `opt` to both solvers. Changing geometry or normalized options invalidates
this cache explicitly. Equal-height, identical grids share one pair of body
factorizations under lateral translation. Different heights require separate
wall blocks. Preparation stores both resistance and mobility factors;
reported setup/factorization costs are shared, not separately incurred by
every subsequent load. Blake caches geometry checks, not FMM trees.

## Geometry and independent diagnostics

`grow_cluster_wall(P,a,delta_body,delta_wall,seed,mode)` returns centres for a
wall at z=0. Translate z for a nonzero `wall_height`. In `"equalHeight"` mode,
all heights are `a+delta_wall`; the lateral configuration depends only on the
seed and body gap. `"threeDimensional"` is a distinct mode, translated to
achieve the minimum requested wall clearance.

`delta_body` is an absolute surface-to-surface distance. Each inserted sphere
has a neighbour at that gap; every pair satisfies the lower bound. Geometry
tolerance is `1e-9*max(a,delta_body)`. The result reports actual minimum gaps
and nearest-neighbour gaps. A singleton has `minBodyGap=NaN`.

`wall_surface_validation(problem,lambda,U)` evaluates a separate surface
quadrature on every sphere. **No special wall-facing points or patches are
added.** Its error is maximum pointwise velocity discrepancy divided by the
global maximum rigid speed; zero-speed cases report absolute error instead.
`linearSystemResidual` in solve details is the physical collocation discrepancy;
`trueGmresResidual` is the independently recomputed square-system residual.

`wall_force_torque_check(problem,lambda)` integrates complete direct Blake
traction. `blake.fields` returns `sigma*n` for outward particle normals; the
diagnostic negates this to obtain the force-on-fluid sign. It sums all sources,
including other particles, and compares integrated force and torque with
physical-source moments. This is a small-reference diagnostic, not a mandatory
step in every solve. Quadrature refinement is needed before trusting it.

A development check at `delta_wall/a=0.1`, `a=mu=1`, requested N=700
(actual N=686, M=801), Rp=0.63 and cutoff 1e-12 used
`U=[1;0.2;0.3;0.1;0.2;0.1]` with direct Blake evaluation. Independent
`[46,56]` surface quadrature gave resistance/mobility boundary errors
`9.315e-3` / `9.306e-3`, despite a relative 2-norm round-trip velocity error
of `3.496e-6`.
Thus this default-resolution narrow-gap case is **unresolved** against the
1e-6 surface target. A successful round trip does not certify surface accuracy.

## Dated experiment

`experiments/okt26_wall_gap_sweep.m` is a script with fixed parameters at the
top: P=10, a=mu=1, seed 20261026, body gap/a=1, N=700, Rp/a=0.63, and eleven
logarithmic wall gaps from 10a to 0.1a. It uses the wall preconditioner, FMM,
GMRES tolerance 1e-9, FMM tolerance 1e-12, and SVD cutoff 1e-12.

The lateral centres and random resistance/mobility data are held fixed over the
sweep. Each wall gap gets one resistance solve and one mobility solve for all P
particles. The script leaves `wallGaps`, the input data, solution arrays,
iteration counts, and residuals in the workspace; it does not save files or
create figures.

```matlab
run('experiments/okt26_wall_gap_sweep.m')
```

The no-argument self-test of `solve_resistance_wall` runs a small resistance
solve and mobility solve, then plots the two-sphere geometry.

## Random unit-torque cluster demo

`demo/demo_random_wall_torque.m` defaults to 100 randomly packed spheres with
unit z-torque on every sphere, a minimum body gap of `a`, and a minimum wall
gap of `0.5a`. It shows the seeded geometry before solving, then visualizes
rigid-body surface speed, independent relative velocity residual, and Blake
traction magnitude. It uses the FMM Blake backend and the free-space one-body
preconditioner; run from the project root with:

```matlab
setpath;
blake.setup;
results = demo_random_wall_torque;
```

For a quick small run, call `demo_random_wall_torque(5)` instead. The returned
structure includes the centres, prescribed torque, rigid velocities, traction
magnitudes, solver details, geometry metadata, and both figure handles.

## Focused validation

```matlab
report = test_okt26_wall_solvers;        % Includes compiled FMM checks.
report = test_okt26_wall_solvers(false); % Direct-only; FMM remains incomplete.
```

Tests cover explicit versus matrix-free square operators, zero loads, cache
invalidation, source moments, actual-height factor reuse, geometry modes,
dimensional scaling, distant-wall drag, small dense comparisons, and an FMM
fixture with 1200 sources. A tiny unsaved fixture also exercises reference
tables, spatial refinements and integrated hidden plotting. It does not run
the production wall sweep. Computational functions have individual no-argument
self-tests. Each self-test prints the quantities it checks and a success
message; geometry/grid and traction self-tests also leave figures open for
inspection. For a batch run, set `groot.DefaultFigureVisible='off'` before
calling the checks. Run MATLAB Code Analyzer on changed files as well.
