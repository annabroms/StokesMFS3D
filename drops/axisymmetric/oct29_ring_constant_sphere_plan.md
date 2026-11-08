# Plan: ring-constant Cartesian Stokeslets on a sphere (MFS note, Section 10, updated PDF `doc/stagnant_cap_mfs-8.pdf`)

Repository: `/home/anna/Dokument/StokesMFS3D` (the `Dokument` copy is the current one).
All new files go in `drops/axisymmetric/`. The note (`doc/stagnant_cap_mfs-7.pdf`) is not touched,
and no existing solver or default behaviour changes.

## Problem and approach

One sphere (radius `a`, centre `c`) in free space with prescribed force `Fg = [0;0;Fz]`
and unknown translation `U`. Stokeslets sit only on a concentric Gauss–Legendre (GL) proxy sphere.
Ring `r` has `n_r` nodes that all carry the same Cartesian vector `a_r` (`lambda_rj = a_r`, no
azimuthal rotation). Two solves on identical geometry and operators:

- `phi = 0`, clean bubble: normal velocity relative to `U` is zero and tangential traction is zero.
  Reference `U = Fg/(4 pi mu a)`.
- `phi = pi`, rigid sphere: full velocity `u = U` with no traction rows.
  Reference `U = Fg/(6 pi mu a)`.

Force is imposed by completion and projection in the updated Section 10.3 (eqs. 74–80), with no
nullspace basis `Z`. Unknowns are auxiliary ring coefficients `alpha in R^{3R}`:

- `Pi_R = I - K_R W'/N`, `a0 = K_R Fg/N`.
- `a = Pi_R alpha + a0`, `U = -beta W' alpha`, `beta = 1/(6 pi mu a)`.
- Boundary equations are written in `alpha`, `A_phi Pi_R alpha + beta V_phi W' alpha = -A_phi a0` (eq. 78), size `3M x 3R`.
- Equivalently the rank-three update `A_phi + (beta V_phi - A_phi K_R/N) W'` (eq. 79). `Pi_R` is never formed densely.
- Velocities, tractions, and all residual rows are the kernel matrices applied to the unknown `alpha` through this map. They are not applied to a separate nullspace coordinate.
- Recovery (eq. 77): `t = sum_r n_r alpha_r`, `a_r = alpha_r - t/N + Fg/N`, `U = -beta t`.
  `a_r` is the physical vector plotted and copied to every source. `alpha_r` is never plotted.

## Existing routines to reuse (none modified)

| Need | Routine |
|---|---|
| GL-in-`cos(theta)` plus uniform-azimuth sphere | `geometry/ellipsoids/setupsurfquad.m` → `setupspherequad.m` with `ellipsoid_param(r,r,r)`. Gives `s.x` (3×N), `s.Nu` (nodes per ring), `s.Nv`, `s.v`. Ordering is fast in azimuth, slow in `v`, so ring id = `repelem(1:Nv, s.Nu)`. |
| Spherical-design collocation | `drops/drop_w_stagnant_cap/stagnant_cap_design.m` |
| Independent validation grid (rotated design, ≥5× collocation, Fibonacci fallback) | `stagnant_cap_validation_grid.m` with `phi=0` (no edge rings) |
| Stokeslet velocity and traction kernels | `images/stokes_mfs_velocity_matrix.m`, `images/stokes_mfs_traction_matrix.m` (basis struct, family `"S"`, interleaved `[x1;y1;z1;...]`) |
| Normal map `Nm`, tangent map `Dm`, velocity `G` (÷`muOut`), traction `Tm` | `drops/bubble_w_stagnant_cap/bubble_boundary_operators.m`. Feed it an all-mobile minimal geometry and options with `n_level=0`; this is checked first. Fallback: call the two kernels directly and copy its reference-axis tangent convention. |
| Field evaluation (velocity, traction, normal velocity, tangential traction) | `bubble_evaluate.m` with `solution = {basis, lambda, U}` |
| Truncated SVD solve | `common/dense/getPseudoFactors.m` |
| Plot conventions | `tiledlayout`, `turbo`, `scatter3` as in the `bubble_plot_*` routines. Those routines are tied to cap and cone structs, so only their style is reused. |

## Files (all in `drops/axisymmetric/`)

Each function header documents purpose, input and output sizes, coefficient ordering, normalization
and coordinate conventions. Each has a no-argument self-test that prints `<name> self-test passed.`
and runs only algebraic or geometric checks, with no physical solve.

1. **`oct29_ring_constant_sphere.m`**: script, run after `setpath`. Sections: PARAMETERS,
   GEOMETRY, ASSEMBLY, SOLVE (two cases), VALIDATION, REPORT, FIGURES.
2. **`ring_gl_sphere.m`**: `[Y, ringId, nr, theta] = ring_gl_sphere(radius, centre, Nu, Nv)`.
   - Wraps `setupsurfquad` and returns proxy points `N×3`, ring membership `N×1`, ring counts
     `R×1`, and ring colatitudes `R×1`.
   - Self-test: ring ids are consecutive, `n_r ≥ 3`, `z = c3 + radius*v`, radius is exact,
     first node of every ring is in the `y=0` meridional plane, and `sum(nr)=N`.
3. **`ring_constant_maps.m`**: `[C, W, KR] = ring_constant_maps(ringId)`.
   - `C = kron(L, I3)` (sparse, `3N×3R`), `W = kron(nr, I3)` (`3R×3`), `KR = kron(ones(R,1), I3)`.
   - Also exposes the implicit projector action `project(alpha) = alpha - KR*(W'*alpha)/N` (eq. 74)
     as a small local function. It is O(R) and the dense `Pi_R` is never built.
   - Self-test: `C'*kron(ones(N,1),I3) = W` (eq. 68), `W'KR = N I3`, `W'*project(x) = 0`,
     `project` is idempotent, and `Pi_R` has rank `3R-3` on a tiny case (dense check only inside the test).
4. **`ring_complete_solve.m`**: `[alpha, a, U, info] = ring_complete_solve(A, V, W, KR, a0, N, beta, rowScale, colScale, tol)`.
   - Builds eq. 79 directly: `M = A + (beta*V - A*KR/N)*W'` (`3M×3R`, no dense projector),
     right-hand side `-A*a0`. Applies `rowScale` and `colScale`.
   - Prints the matrix dimensions before the SVD and calls `getPseudoFactors(·, tol, 0)`.
   - Recovers (eqs. 75, 77): `t = W'*alpha`, `a = alpha - KR*t/N + a0`, `U = -beta*t`. It returns
     `alpha`, `a` (`3R×1`, physical per-node vectors), `U`,
     `info = {sigmaMax, threshold = tol*sigmaMax, rank, size}`.
   - Self-test: a small random system built from known `(a, U)` satisfying eq. 75 is recovered exactly.
     Eq. 79 is checked against the explicit `A*Pi + beta*V*W'` form, `W'a = Fg` holds after
     deliberately truncating the SVD, and the rank and threshold are reported.
5. **`ring_plot_residuals.m`**: `fig = ring_plot_residuals(points, residuals, labels, visibility)`.
   One `1×3` tiled figure of `scatter3` maps at validation nodes, with `log10` residual
   and the normalization named in each colorbar label.
6. **`ring_plot_sources.m`**: `ring_plot_sources(ax, Yrep, a, centre, radius, proxyRadius, arrowScale, titleText)`.
   - Draws physical and proxy meridional sections at `y=0`.
   - Draws one `quiver3` arrow per ring at its `j=0` node (colatitude `theta_r`, `psi=0`) with
     `autoscale off`, oblique view `view(35,22)` so `a_y` is visible, and a text line with the arrow scale.
   - Self-test: invisible figure, `R` arrows present and axes limits finite.

No wrapper layers. `Nu` and `Nv` are the only GL controls.

## Parameters at the top of the experiment

`a=1, mu=1, c=[0 0 0], Fg=[0;0;-1]`, proxy radius `ap=0.6*a`, `Nu=24, Nv=12` (counts rounded
up to even by `setupspherequad`), collocation request `Mreq` (spherical design, e.g. 600),
validation multiplier 5, relative SVD tolerance `1e-10`, quiver fraction 0.3. All values are exposed
and documented.

## Dimensions (symbols; numbers printed at run time)

`R=Nv` rings, `N=sum(n_r)` sources, `M` collocation points, `Mv≥5M` validation points.

| Quantity | Size |
|---|---|
| `G`, `Tm` (nodal, from existing kernels) | `3M×3N` |
| `C` / `GR=G*C`, `TR=Tm*C` | `3N×3R` / `3M×3R` |
| `W`, `KR` | `3R×3`, `3R×3` |
| `Pi_R` (applied implicitly, rank `3R-3`) | `3R×3R` (never formed) |
| Bubble `A=[Nm*GR; Dm*TR]`, `V=[Nm*KM; 0]` | `3M×3R`, `3M×3` (rows `M` normal + `2M` traction) |
| Rigid `A=GR`, `V=KM` (`KM=kron(ones(M,1),I3)`) | `3M×3R`, `3M×3` |
| Completed system `A Pi_R + beta V W'` (eq. 78) | `3M×3R`, unknown `alpha`; effective rank ≤ `3R` |

The scalar unknowns are `3R` (`alpha`). Nodal `3N` unknowns are never solved for, only used to evaluate kernels.

## Completion, scaling, reconstruction

- Velocities and tractions are evaluated as kernels applied to the unknown `alpha` via `a = Pi_R alpha + a0`
  and `U = -beta W' alpha`. The unknown subspace is `range(Pi_R)` of dimension `3R-3`; the three
  force-carrying components of `alpha` encode translation (eq. 75).
- `a0 = K_R Fg/N` gives every Stokeslet strength `Fg/N`, so `W'a0 = Fg`, and `W'Pi_R = 0`.
  Force balance therefore holds to round-off even after SVD truncation.
- Physical reconstruction is `a = Pi_R alpha + a0`, `lambda = C*a` (`3N`, per-node strengths).
  The plots and force checks use this, never `alpha` or `a0` alone.
- Row scaling: velocity rows divided by that case's `U_ref`, and bubble traction rows divided by
  `mu*U_bubble/a`. This makes all rows dimensionless and unit scale.
- Column scaling: `alpha = |Fg| * alpha_hat` with `|Fg|` the force-magnitude unit, so the unknown is dimensionless.
  `beta = 1/(6 pi mu a)` is only a conversion scale, not the mobility, and is the same in both cases. `U` comes from `-beta W' alpha`.
  The relative truncation threshold `tol*sigma_max` applies to this scaled matrix.

## Validation and report

For each case, evaluate on the independent `Mv` grid through `bubble_evaluate(...)`, which uses the
full Cartesian kernels with no azimuthal-symmetry assumption:

- Computed `U` and relative error `|U - U_exact|/|U_exact|` (vector norm, so transverse drift counts).
- Maximum normalized residuals:
  - bubble `|(u-U)·n|/|U_bubble|` and `|(I-nn')t|/(mu |U_bubble|/a)`;
  - rigid `|u-U|/|U_rigid|`.
- Relative force residual `|sum n_r a_r - Fg|/|Fg|` (via `W'a`) and `|sum_j lambda_j - Fg|/|Fg|`.
- Hydrodynamic torque `|sum (y-c)×lambda|/(|Fg| a)` as a diagnostic (note eq. 86).
- Printed before and after SVD: `R, N, M, Mv`, unknowns `3R`, matrix size, `sigma_max`, threshold, rank.

The report states honestly that this is a restricted approximation space. A good `U` with large
residual maps is a negative result, not something to tune away.

## Figures

1. `ring_plot_residuals`: bubble normal-velocity, bubble tangential-traction and rigid
   full-velocity residual maps on the validation nodes, each labelled with its normalization.
2. `ring_plot_sources`, two panels (bubble, rigid) in one `tiledlayout`. They show the physical shared
   vector `a_r` per ring, one location per ring at `psi=0`. Both use one shared arrow scale
   `0.3*ap/max_r|a_r|` over both cases, stated on the figure as per-node force units, with the
   physical sphere as a solid section and the proxy sphere as a dashed one.

## Todos (also tracked in SQL)

0. `save-plan`: first step after approval, save this full plan as `drops/axisymmetric/oct29_ring_constant_sphere_plan.md` (separate from the note), then implement.

1. `ring-gl-sphere`: write `ring_gl_sphere.m` plus self-test.
2. `ring-constant-maps`: write `ring_constant_maps.m` (`C`, `W`, `KR`, projector action) plus self-test.
3. `ring-complete-solve`: write `ring_complete_solve.m` (eq. 79 rank-three update, recovery eq. 77) plus self-test.
4. `ring-plot-residuals`: write `ring_plot_residuals.m` plus self-test.
5. `ring-plot-sources`: write `ring_plot_sources.m` plus self-test.
6. `ring-experiment`: write `oct29_ring_constant_sphere.m`, which first confirms the
   `bubble_boundary_operators` adapter and the `stagnant_cap_validation_grid(phi=0)` call, then
   assembles `G`, `Tm`, reduces by `C`, runs the two solves, validates and plots.
7. `ring-verify`: run the five self-tests and the experiment in MATLAB (R2025a is installed at
   `/usr/local/MATLAB/R2025a/bin/matlab`). Check that the printed dimensions match the table,
   force residual ≈ 1e-13, and that figures render. Confirm that `git status` shows only
   `drops/axisymmetric/` additions.

## Notes and risks

- The exact bubble field has radial Cartesian components that vary with azimuth. The tied-Cartesian
  space cannot represent them exactly, so large bubble residuals are possible. The experiment reports
  them and does not hide them.
- Per-node strengths are plotted (the note's convention 63). A ring's total force is `n_r a_r`.
  The figure labels say so.
- Collocation uses a spherical design while the proxy uses GL rings, so the two grids are not aligned.
- Poles are avoided automatically because GL nodes exclude `v=±1`.
- Eq. 81 (explicit `Z = Zs ⊗ I3` via a Householder reflection) is not used. The projector form (eqs. 74–79) is the single formulation.
