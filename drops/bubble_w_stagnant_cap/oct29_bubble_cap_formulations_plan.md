# Exterior-only stagnant-cap bubble: explicit force and translation plan

Status: updated before implementation, 29 October 2026. The force-completion
revision below supersedes the weighted-force-row approach for the primary
comparison, following the user's request. The mathematical document is
unchanged. Preserve the existing working-tree edits and saved experiment
results.

## Current revision: force completion for the Cartesian cap formulation

### What changes and what is reused

The current saved experiment demonstrates why force balance must be imposed
algebraically. Its baseline explicit-force least-squares solve has relative
force residual `0.925235`, with `Uz=-0.02226496` against theoretical
`Uz=-0.30780528`. The previous solver has `Uz=-0.29780218` and force residual
`7.39e-13`. The smaller boundary residuals and nearly uniform surface maps of
the under-forced solution do not establish better accuracy.

Reuse [bubble_system](./bubble_system.m) and [bubble_solve](./bubble_solve.m),
adding an optional `capMode="cartesian"` argument; the default remains
`"split"`. Both modes share the same force projection, completion, source
basis, translation recovery, SVD, and field evaluation. No new solver helper
or duplicated kernel/geometry routine is needed. Retain
[bubble_explicit_system](./bubble_explicit_system.m) and
[bubble_explicit_solve](./bubble_explicit_solve.m) as the existing direct-force,
weighted least-squares alternative, clearly labelled as such.

### Force-free physical sources and the known completion flow

For `N` Stokeslet sites, including cone sites, define

\[
H=K_N^T\in\mathbb R^{3\times3N},\qquad
C=H^T(HH^T)^{-1}\in\mathbb R^{3N\times3},\qquad
P_F=I-CH,\qquad \lambda_0=CF_g.
\]

Use the existing algebraic parametrization

\[
\lambda_{free}=P_F\alpha,\qquad
\lambda=\lambda_{free}+\lambda_0,\qquad
U=-\beta H\alpha,\qquad \beta=\texttt{completionScale}.
\]

Then `H*lambda_free=0`, `H*lambda0=Fg` and `FH=-H*lambda=-Fg`, up to
floating-point cancellation, independently of the truncated-SVD fit. All
unknown **physical** source contributions are force-free; the known flow
`G*lambda0` carries the prescribed force. The auxiliary `alpha` is not itself
the physical force-free vector: its three force-map components encode the
independent bubble translation. This is exactly the existing solver's
parametrization and avoids adding a redundant independent translation vector.

Apply projection by low-rank products, e.g. `G-(G*C)*H`; never form the dense
`3N x 3N` projector. The completion contributes to every boundary right-hand
side through velocity or traction evaluated by the existing kernels. There
are no fitted force rows and no force penalty weight in this revised solve.
Continue reporting the actual physical force residual and raw defect.

### Cartesian boundary matrix and dimensions

Keep interleaved nodal Cartesian components. With `M=Mc+Mm`, use

\[
\begin{bmatrix}
G_cP_F+\beta K_{M,c}H\\
\mathcal N_mG_mP_F+\beta B_mH\\
D_mT_mP_F
\end{bmatrix}\alpha
=-\begin{bmatrix}
G_c\lambda_0\\
\mathcal N_mG_m\lambda_0\\
D_mT_m\lambda_0
\end{bmatrix},\qquad B_m=\mathcal N_mK_{M,m}.
\]

The row counts are `3Mc`, `Mm`, `2Mm`, respectively, so the revised matrix
has size `3M x 3N`, matching the previous solver. `Gc/Gm` map source forces
to physical velocity (`1/muOut` included); `Tm` maps them to physical traction.
`KM,c/KM,m` repeat translation, `Nm` extracts the normal component, and `Dm`
is the existing two-row orthonormal mobile tangent representation. Impose
full Cartesian velocity on the cap, with no separate cap normal rows.

The old split-cap rows are an orthogonal row transformation/permutation of
these Cartesian rows when each cap component has the same point weight.
Consequently the two formulations have the same least-squares objective and
singular values in exact arithmetic under matched scaling. Compare physical
fields and translation; do not require redundant coefficient vectors to agree.

### Scaling and diagnostics

The comparison experiment uses only one configured `svdTolerance` (default
`1e-8`), shared by both modes and every resolution. Do not run a cutoff sweep.
Its agreement budget uses measured resolution changes only and explicitly
states that cutoff sensitivity has not been assessed. Multiple-cutoff checks
in the small solver self-tests remain separate from the experiment.

Reuse existing physical scales `U*=norm(Fg)/zetaBubble` and
`traction*=muOut*U*/a` for both modes. Divide velocity rows by `U*` and traction
rows by `traction*`; use the same `beta` and relative SVD cutoff. The primary
experiment keeps unit condition weights with area, column and condition-count
scaling off. Permit area and column scaling through existing code; Cartesian
cap rows require equal **effective** normal/tangential component weights at
each cap target. Reject incompatible anisotropic weights rather than silently
changing the objective.

Print the mode, three boundary row-block counts, `3M x 3N` dimensions and
memory preflight before SVD. Retain relative/absolute threshold, rank, and
physical force-defect reporting. Verify force balance for different cutoffs,
not just for a single retained rank.

### Implementation sequence (reuse existing functions)

1. **Update this plan before code edits.** This revision records the agreed
   force-completion design separately from the mathematical document.
2. **Assembly:** extend `bubble_system` with an optional fourth `capMode`
   argument. Default `"split"` preserves the existing matrix and RHS exactly.
   `"cartesian"` reuses its projected kernels and known completion, changing
   only the cap/mobile boundary row grouping. Return explicit row metadata.
3. **Solve/scaling:** extend `bubble_solve` with an optional fifth `capMode`
   argument (existing fourth `pointWeights` argument unchanged). Branch only
   row scales/weights and row labels. Reuse the same pseudoinverse and physical
   recovery. Return the mode, `lambda_free`, `lambda0`, and raw force defect so
   the prescribed-force and force-free contributions can be audited.
4. **Tests:** extend the existing end-of-file no-argument self-tests rather
   than add separate one-use helpers. Check dimensions, local cap-row
   equivalence, completion RHS, force sign and zero force of `lambda_free`,
   endpoint empty blocks, cone nodes, and matching physical fields/translation
   under multiple cutoffs. Exercise supported area/column flags and reject
   unequal effective cap-component weights. Recheck the saved legacy baseline.
5. **Experiment:** update
   [oct29_bubble_cap_formulations](./experiments/oct29_bubble_cap_formulations.m)
   to compare `"split"` and `"cartesian"` modes of `bubble_solve`. Keep all
   comparison calculations and plotting in the script. Both modes use the
   same geometry, physical scales, independent grids and solver settings.
   Remove force-weight controls from this primary comparison; both matrices
   are `3M x 3N`. Add force-free/completion diagnostics and unambiguous method
   labels to exported metadata. Keep the weighted direct solver available
   for separate diagnostic use.
6. **Plots and outputs:** add a dedicated frame/plot explanation, recovered
   `Uz/U*` and force residual to the surface figures; visibly mark the cap edge.
   Preserve the signed polar component and total relative speed, explicitly
   noting that total speed includes normal error. The existing separate
   normal-residual map diagnoses that error. Save revised outputs under a
   distinct force-completion prefix, preserving the existing direct-force
   results for comparison. Do not silently relabel or overwrite those results.
7. **Validation:** run Code Analyzer and small no-argument checks. Updating the
   script does not execute the production experiment or regenerate its figures.

### Reading the bubble-relative surface-velocity figure

Every panel uses `v=u_out-U` with each solution's own recovered translation;
`U*` is the same nonzero theoretical reference speed for both solutions.
The cap is `theta<phi`, the upper spherical patch for the configured polar
axis. On it, `v=0` is imposed. On the mobile region, only `n dot v=0` and
zero tangential traction are imposed: nonzero tangential velocity is expected.

The upper row shows signed `v dot e_theta/U*`, positive towards increasing
polar angle (from the north pole towards the south pole). Negative color
means motion in the opposite direction; it is not a residual magnitude.
The lower row shows `norm(v)/U*`, which contains polar, azimuthal and normal
components. Only for axisymmetric, no-penetration flow does it reduce to
`abs(v dot e_theta)/U*`. Nonzero cap values measure boundary error; mobile
values mostly represent allowed slip. Normal leakage is checked separately.

In the saved weighted-force experiment the explicit solution moves much too
slowly, so its small relative-velocity plot primarily reflects under-forcing.
Use the recovered translation, theoretical speed error and force balance
alongside this figure. For the revised paired projected formulations, expect
closely matching fields under equivalent row weighting, subject to finite
precision and truncation sensitivity.

## Historical direct-force-row design

The sections below record the original weighted-force-row implementation.
They are retained as the specification for the optional direct-force solver;
the force-completion revision above governs the updated comparison experiment.

## 1. Scope and invariants

Add an exterior-only solver for independent physical Stokeslet strengths and
translation, `(lambda,U)`. In the laboratory frame the exterior fluid is at
rest at infinity. Impose

\[
u_{out}=U\quad\hbox{on the cap},\qquad
n\cdot(u_{out}-U)=0,\quad \Pi\sigma_{out}n=0
\quad\hbox{on the mobile surface},\qquad \Pi=I-nn^T.
\]

Use three Cartesian velocity rows per cap point, with **no additional cap
normal rows**. Retain the existing two orthonormal tangent rows for mobile
traction. Include three force rows. Do not project the new coefficients,
insert a known-force completion, eliminate translation, or repair force balance
after solving. There is no interior field, background uniform flow, angular
velocity unknown, torque equation, or normal-traction condition.

Keep [bubble_solve](./bubble_solve.m) and its default behavior available as the
previous split-cap, projection-and-completion formulation. Its options, physical
recovery, default cutoff, and existing experiments remain compatible. Shared
code extraction must preserve its matrix, right-hand side, scaling and results.

The new path and comparison use only unit-scale `S` blocks:
`refinementFamilies="S"`, `derivative_direction="none"`. Optional cone Stokeslets
remain supported through the existing geometry options. Assert
`coefficientCount == 3*sourceCount`; here `N` means actual Stokeslet sites,
including base and cone sites. Preserve existing placement, clustering,
collocation patches, azimuth counts, and refinement rules. Other families remain
available in the previous solver; extending the new formulation to them is
outside this comparison.

## 2. Inspected implementation and reuse

| Existing routine | Reuse and relevant observed behavior |
| --- | --- |
| [bubble_options](./bubble_options.m) | Completes physical/geometry controls; default `svdTolerance=1e-12`, `conditionWeights=[1 1 1]`, with area, column and condition-count scaling off. Keep these defaults. |
| [bubble_geometry](./bubble_geometry.m), [bubble_basis](./bubble_basis.m) | Build interior base/cone sites, target masks, normals, counts and basis ordering. Geometry delegates to [stagnant_cap_geometry](./stagnant_cap_geometry.m), [stagnant_cap_refined_geometry](./stagnant_cap_refined_geometry.m), and [stagnant_cap_cluster_offsets](./stagnant_cap_cluster_offsets.m). Build once per comparison pair. |
| [bubble_source_refinement](./bubble_source_refinement.m) | Existing optional source-depth collocation patches; keep its rules and pass the resulting identical target cloud to both methods. |
| [bubble_system](./bubble_system.m) | Existing `M` normal, `2Mc` cap tangential, `2Mm` mobile traction rows. Extract its unprojected kernels, sparse normal/tangent maps and index lists for common use. Preserve its tangent-frame construction exactly. |
| [stokes_mfs_velocity_matrix](../images/stokes_mfs_velocity_matrix.m) | Laboratory velocity kernel; divide by `muOut` for physical strengths. Includes the Stokeslet `1/(8*pi)` factor. |
| [stokes_mfs_traction_matrix](../images/stokes_mfs_traction_matrix.m) | Physical traction on the drop-to-exterior normal; no viscosity multiplier for fixed physical Stokeslet strengths. Its existing self-test verifies integrated traction has force `-lambda`. |
| [getKmat](../common/dense/getKmat.m) | Reuse its first three columns for source translation/force maps; last three columns permit an optional source-torque diagnostic. No rotation unknown is needed. |
| [image_mobility_maps](../images/image_mobility_maps.m) | Retain in the previous solver. It computes a completion even when only its first output is requested, so the new solver will use `getKmat` instead. |
| [getPseudoFactors](../common/dense/getPseudoFactors.m), [bubble_svd_memory](./bubble_svd_memory.m) | Economy truncated SVD and dense-memory preflight. Retained singular values satisfy `sigma > tau*sigmaMax`; the third output is a vector of retained values. |
| [bubble_validate](./bubble_validate.m), [bubble_validation_grid](./bubble_validation_grid.m), [bubble_check_validation](./bubble_check_validation.m) | Shared physical residual definitions, shifted independent nodes, edge refinement and second-grid checks. Preserve the existing four-entry `report.maxima` layout. |
| [stagnant_cap_validation_grid](./stagnant_cap_validation_grid.m) | Rotates the base design by fixed nonzero angles and adds azimuthally shifted edge rings; falls back to a Fibonacci grid for large requests. Explicitly check independence from collocation. |
| [stagnant_cap_exact_drag](./stagnant_cap_exact_drag.m) | Bubble theory is the `muIn=0` argument of this analytic function; do not set an interior solver viscosity to zero. |
| [bubble_plot_geometry](./bubble_plot_geometry.m), [bubble_plot_residuals](./bubble_plot_residuals.m) | Reuse geometry views and residual-map conventions. Add optional axes/limits/condition selection to the residual plotter while preserving existing calls. |
| [plot_stagnant_cap_velocity](./plot_stagnant_cap_velocity.m) | Reuse its polar tangent, relative-velocity and edge-profile conventions. Its current API is for two-phase options and creates its own grids/figures, so it cannot directly supply the paired common-grid plots. |
| [oct24_stagnant_cap_bubble_residuals](./experiments/oct24_stagnant_cap_bubble_residuals.m) | Reuse its physical setup, Stokeslet configuration, preflight, shared validation and results-export pattern. |

The generic [drop_plane_velocity](./drop_plane_velocity.m) selects interior and
exterior formulas for two-phase drops. Use the shared exterior evaluator below
for the masked bubble plane instead of introducing interior data or changing
this routine. No kernel or geometry implementation will be duplicated.

## 3. Explicit block system

All nodal Cartesian arrays are interleaved:
`[x1;y1;z1;x2;y2;z2;...]`. Basis blocks are concatenated in the order returned
by `bubble_basis`: base Stokeslets, then cone Stokeslets. Within each surface
region, preserve the original target-index order. Unknowns are

\[
x=[\lambda_{x1},\lambda_{y1},\lambda_{z1},\ldots,
\lambda_{xN},\lambda_{yN},\lambda_{zN},U_x,U_y,U_z]^T.
\]

Let `M=Mc+Mm`. With the existing exterior velocity and traction kernels, form

\[
\underbrace{\begin{bmatrix}
G_c&-K_{M,c}\\
\mathcal N_mG_m&-B_m\\
D_mT_m&0_{2M_m\times3}\\
K_N^T&0_{3\times3}
\end{bmatrix}}_{A\;\in\;\mathbb R^{(3M+3)\times(3N+3)}}
\begin{bmatrix}\lambda\\U\end{bmatrix}
=
\underbrace{\begin{bmatrix}
0_{3M_c}\\0_{M_m}\\0_{2M_m}\\F_g
\end{bmatrix}}_b,
\qquad B_m=\mathcal N_mK_{M,m}.
\]

| Block | Dimensions | Meaning |
| --- | --- | --- |
| `Gc`, `Gm` | `3Mc x 3N`, `3Mm x 3N` | Exterior physical velocity at cap/mobile targets; `G=stokes_mfs_velocity_matrix(basis,points)/muOut`. |
| `KM,c`, `KM,m` | `3Mc x 3`, `3Mm x 3` | Repeat a Cartesian translation: `kron(ones(Mregion,1),eye(3))`. |
| `Nmap_m` (`mathcal N_m`) | `Mm x 3Mm` | One local `n_i'` row at each mobile point. |
| `Bm` | `Mm x 3` | Normal component of repeated translation; row `i` is `n_i'`. |
| `Tm` | `3Mm x 3N` | Exterior `sigma*n` at mobile targets with outward unit normals. |
| `Dm` | `2Mm x 3Mm` | Existing mobile `t1'`, `t2'` rows, interleaved per target; local `D_i'*D_i=I-n_i*n_i'`. |
| `KN` | `3N x 3` | Repeat translation at every Stokeslet source; obtain from the first three columns of `getKmat(sourcePoints,centre)`. |
| `KN'` | `3 x 3N` | Sum physical Cartesian source forces, including all cone Stokeslets. No source-area quadrature factors. |

Row blocks in order: cap Cartesian `3Mc`, mobile normal `Mm`, mobile tangent
traction `2Mm`, force `3`. Assert their sum is `3M+3`. Keep the zero translation
columns in the traction and force blocks. Endpoint tests allow empty cap or
mobile blocks without adding substitute equations.

The force sign is

\[
F_H=-K_N^T\lambda,\qquad F_H+F_g=0
\quad\Longleftrightarrow\quad K_N^T\lambda=F_g.
\]

Return physical `lambda`, `U`, `FH`, `Fg`, basis/grid metadata, scales,
matrix dimensions, row/column scalings, singular cutoff and retained rank in
the same solution-field conventions as the previous solver. Recover `U`
directly from the last three unknowns. Optional torque remains a diagnostic.

For comparison, the previous Stokeslet-only system is `3M x 3N` for `alpha`,
with `K=KN'`, `C=K'/(K*K')`, `lambda0=C*Fg`,
`lambda=alpha-C*(K*alpha)+lambda0`, and `U=-beta*K*alpha`.
Its rows are ordered `[all normals; cap tangents; mobile tangent traction]`.
Document `beta=completionScale`, default `1/(6*pi*muOut*a)`, in the output;
it is absent from the new solve.

## 4. Scaling, weighting and SVD policy

For nonzero axial `Fg`, use the same physical references for both solutions:

\[
\zeta_b=\texttt{stagnant\_cap\_exact\_drag}(\phi,0,\mu_{out},a),\quad
U_{th}=F_g/\zeta_b,\quad U_*=\|U_{th}\|>0,\quad
t_*=\mu_{out}U_*/a>0,\quad F_*=\|F_g\|>0.
\]

These agree with the previous solver's axial `Uref` and `tractionScale`.
Normalize velocity rows by `1/U*`, traction rows by `1/t*`, and new force
rows by `wF/F*`, with explicit dimensionless `wF=1` in the primary experiment.
The scaled force right-hand side is `wF*Fg/F*`. Force rows have no surface-area
weight and no collocation-count divisor. Both methods minimize the resulting
scaled residual in least squares; report unweighted physical residuals as well
as weighted fit norms.

Primary comparison settings: `conditionWeights=[1 1 1]`,
`useConditionCountScaling=false`, `pointWeights=ones(M,3)`,
`useAreaWeights=false`, and `useColumnScaling=false`, explicitly set in the
script. This preserves the inspected experiment's scaling choices and keeps
cap weighting isotropic: all three components at a cap point have one weight.
For `r=u-U` and `Q_i=[n_i';t1_i';t2_i']`,

\[
\|Q_i r\|_2^2=r_n^2+r_{t1}^2+r_{t2}^2=\|r\|_2^2.
\]

Consequently replacing split cap rows by Cartesian rows preserves their
least-squares contribution before the formulation's other changes. This
identity applies with a common pointwise scalar area weight too. Unequal
normal/tangent weights do not preserve it: in the existing options the three
condition weights refer to normal velocity, cap tangential velocity and mobile
traction. Condition-count scaling divides these by `sqrt(M)`, `sqrt(Mc)` and
`sqrt(Mm)`, respectively, so it generally introduces anisotropic cap weights.
Reject such settings in this comparison rather than silently transferring
them to Cartesian components. General anisotropic weighting would require
the Cartesian block weight `Q_i' * diag(wn,wt,wt) * Q_i`; it is not part of
the primary implementation.

Retain area/column options in the legacy path. The new solver initially accepts
the isotropic, unweighted comparison settings, and fails clearly for unsupported
weighting flags rather than ignoring them. New-solver-only settings, e.g.
`forceWeight`, live in a separate completed settings structure, without
changing `bubble_options` defaults or its allowed option names. If area or
column scaling is later enabled for the new solver, first extract/reuse the
legacy area/scaling code and add equivalence tests. Column scaling, when used,
must include translation columns and be undone on all unknowns.

Use `getPseudoFactors(A_scaled,tau,0)` with `tau=1e-8` for the primary script,
matching the inspected current residual experiment. Keep the library default
`1e-12` unchanged. Before **every** SVD, including tests and sensitivity solves,
print actual matrix dimensions, block counts, dense-memory preflight and `tau`.
Afterward print `sigmaMax`, absolute threshold `tau*sigmaMax`, and retained
rank over `min(nrows,ncols)`. With at least one retained value, the first
returned retained singular value gives `sigmaMax`; detect an empty retained
spectrum explicitly. Do not compute a second SVD merely to obtain the threshold.

Including force rows in overdetermined least squares does **not** make force
balance exact. Always print and save the three-component force defect and
`norm(KN'*lambda-Fg)/F*`, separately from the boundary fit. The previous solver
enforces it algebraically up to floating-point cancellation; the new solver
may trade it against boundary residuals. Neither a small cutoff nor the
presence of force rows is an accuracy certificate. Do not adjust force weight
silently or substitute a constrained solve. If force error is material,
report it and its sensitivity explicitly.

## 5. Proposed files and compatibility changes

The files below implement the reviewed design. Comparison calculations and
plotting remain directly in the experiment script, as requested during review.

| Proposed file | Responsibility and reusable inputs |
| --- | --- |
| `drops/bubble_boundary_operators.m` | Common unprojected kernel assembly, source basis, target subsets, existing tangent frames and sparse normal/tangent maps. Extract from `bubble_system`, preserving ordering. No force projection/completion. |
| `drops/bubble_explicit_system.m` | Assemble the four explicit blocks and metadata using those operators and `getKmat`; assert Stokeslet-only basis and all block dimensions. |
| `drops/bubble_explicit_solve.m` | Complete new settings once; validate allowed weights/flags; apply dimensional row scales, memory guard and SVD; recover independent physical strengths/translation; report force defect and rank. Suggested API `(opt,Fg,geometry,settings)`. |
| `drops/bubble_evaluate.m` | Shared batched exterior velocity and optional traction evaluation from `solution.basis`, physical `lambda`, and `muOut`; optionally project supplied normals. Accept arbitrary exterior/surface points for validation, profiles and the meridional grid. |
| `drops/experiments/oct29_bubble_cap_formulations.m` | Requested script: configuration, one shared geometry, two solves, shared independent validation, per-method metrics and cross-field differences, agreement diagnostics, tables, six visualization groups, and exports. Compute comparisons and assemble plots directly in clearly separated script sections, reusing the shared evaluator and existing geometry/residual plotters. Keep evaluated arrays, figure handles and plot metadata available in the script workspace. No automatic execution from tests. |

Focused compatibility edits:

- [bubble_system](./bubble_system.m): call the extracted boundary operators;
  retain its existing projection/completion assembly and physical recovery.
- [bubble_validate](./bubble_validate.m): call the shared batched evaluator;
  add full cap velocity norm and optional cached physical fields. Preserve
  existing four `maxima` entries, units and default memory behavior. The full
  cap norm is evaluated pointwise, not reconstructed from separate maxima.
- [bubble_plot_residuals](./bubble_plot_residuals.m): optional supplied axes,
  shared limits and condition selection, including full cap velocity. Existing
  calls retain their original four-panel display and color behavior.
- [bubble_solve](./bubble_solve.m): add returned/printed absolute SVD threshold
  using its already returned singular values, without changing the solve.
- Add a short link from [README_stagnant_cap](./README_stagnant_cap.md) once
  implementation exists. Leave the mathematical document and the original
  implementation plan unchanged.

Every new public function dispatches `nargin==0` to a small deterministic
`self_test` at the end of its file. Avoid trivial helper functions; any new
substantive local helper also needs a meaningful no-argument test path and a
one-line purpose comment. Self-test routines themselves are the no-argument
checks. Comparison calculations and plotting belong directly in the October 29
script because they serve this experiment. Keep them in readable sections
without introducing separate comparison/plotting functions or local wrappers
solely for organization. Reusable assembly, solving and field evaluation remain
in the tested functions above.

## 6. Experiment construction and independent accuracy

Start from the current Stokeslet-only configuration in the inspected
[residual experiment](./experiments/oct24_stagnant_cap_bubble_residuals.m):
`a=muOut=1`, centre zero, `phi=pi/4`, `sourceRequest=700`,
`innerProxyRatio=0.63`, `collocationRatio=2`, five optional cone levels,
`n_az_source=64`, `n_az_colloc=80`, and ten collocation rings per region.
Retain source start distance `0.35a`, tapered source clustering with sigma 3,
geometric collocation ratio 0.5 and extent 0.5 radians, and its existing
source-depth patch settings. Copy and save the complete normalized options,
including inactive controls; actual catalogue-rounded counts are authoritative.
Expose the existing cone switch through `n_level`: zero gives the uniform
proxy and a positive value uses the existing cone-node rules. Build a fresh
paired geometry for each chosen configuration. Do not change the refinement
law when switching formulations. For `n_level=0`, set
`useSourceRefinement=false` as required by `bubble_options`; cone footprints
do not exist in that configuration.

Compute `Fg=(4*pi*a^3/3)*densityDifference*gravity(:)`, with the inspected
`densityDifference=1`, `gravity=[0 0 -1]`. Pass the **same geometry object**
to `bubble_solve` and `bubble_explicit_solve`. Assert exact equality of returned
source coordinates, collocation coordinates, normals, region masks, basis
ordering and counts. Do not call adaptive collocation separately per solver.
Complete and record all settings before either solve.

Build one primary validation grid with at least `5*M` points and a finer
audit with at least `10*M` points using existing validation builders. Use
rotated/shifted base points and finer, staggered edge rings on both sides.
Apply the existing finer-offset recipe in `bubble_check_validation`, then
explicitly verify positive edge distances and no collocation coincidences
within `1e-11*a`. If coincidences occur, deterministically shift the validation
azimuths, recheck, and record the shift; never move collocation or source nodes.
Remove any exact-edge validation samples and retain/replenish the requested
off-edge base density. Keep the density requirement after any exclusions.
When source-depth patches
are enabled, audit their independence too rather than assuming all appended
points are shifted.

Define edge distance consistently as surface arc distance
`s=a*abs(theta-phi)`. Print actual minimum `s` and `s/a` for cap/mobile and
overall on each grid, and the closest source and collocation distances.
Exclude the exact edge throughout; do not replace computed edge traces by
their theoretical limit. Report primary and fine-grid maxima and changes;
flag changes above the existing 10% criterion as unresolved sampling. Evaluate
the following metrics for both solutions on the same larger grids:

\[
e_U=\|U-U_{th}\|/U_*,\quad
r_c=\max_{cap}\|u-U\|/U_*,\quad
r_n=\max_{mobile}|n\cdot(u-U)|/U_*,
\]
\[
r_t=\max_{mobile}\|\Pi\sigma n\|/t_*,\qquad
r_F=\|F_H+F_g\|/F_*.
\]

Also report signed `Uz`, all three translation components, transverse drift,
and legacy cap normal/tangential diagnostics. Counts include base/cone/total
sources, physical-strength unknowns, total unknowns, base/ring/patch target
counts, `Mc`, `Mm`, validation counts, matrix dimensions and retained rank.
An absent endpoint condition is `NaN`/not applicable, never an invented zero.

Cross-method quantities use common scales, never either recovered speed as
denominator:

\[
\Delta_U=\|U_{old}-U_{new}\|/U_*,\quad
\Delta_u=\max_{surface}\|u_{old}-u_{new}\|/U_*,\quad
\Delta_t=\max_{surface}\|\Pi(t_{old}-t_{new})\|/t_*.
\]

These field differences use laboratory velocities and common normals. Report
tangential traction difference separately on the cap, mobile surface and
whole surface: cap shear is an unconstrained reaction. Additionally record
the relative-velocity difference
`((u_old-U_old)-(u_new-U_new))/U*`, clearly distinguished from the lab-frame
difference. Use the same definitions on the exterior meridional grid. Source
coefficients are representation-dependent and may be redundant: do not use
their agreement as an acceptance criterion.

## 7. Required plots

Use one set of evaluated arrays for each common grid. Corresponding panels
share coordinates, camera, axes, color map, limits and normalization. Compute
limits from the union of paired values. Signed components have symmetric
limits; nonnegative norms have a zero lower bound. Residual maps may use
`log10(max(residual,1e-16))`, with the same floor and per-condition paired
limits; tables retain unfloored values. Mask inapplicable regions rather than
coloring them as zero. Avoid surface interpolation across the cap edge.

1. **Surface velocity:** side-by-side maps of `(u-U)/U*` using each recovered
   `U`. Show signed meridional component and vector magnitude, with optional
   vector arrows; pairwise common limits for each scalar.
2. **Surface difference:** map `norm(u_old-u_new)/U*` on the shared independent
   surface points. Label it as a laboratory-frame difference; optionally add
   the relative-velocity difference in a separately labelled panel.
3. **Boundary residuals:** paired maps for full cap velocity on the cap,
   mobile normal velocity on the mobile region, and mobile tangential
   traction on the mobile region. Retain paired legacy split-cap normal and
   tangential diagnostic maps as supplementary panels. Thus every old and
   new imposed condition can be inspected, with common limits for matches.
4. **Polar profiles:** evaluate at identical off-collocation polar angles,
   with logarithmic edge offsets on both sides and fixed common azimuths.
   Plot signed `(u-U) dot e_theta/U*` for both methods on the same axes and
   their difference, with full-angle and close edge views. Mark `phi`, break
   lines across the excluded edge, and state the nearest sample. Additional
   azimuth traces or an envelope expose departures from axial symmetry.
5. **Exterior meridional field:** one fixed plane through the bubble centre,
   e.g. `y=centre(2)`, box `[-3a,3a]` in relative x/z coordinates, with common
   near-surface refinement. Set interior values to `NaN` and do not evaluate
   interior sources as a physical interior flow. Plot paired laboratory
   velocity magnitude/components and arrows with shared limits and arrow
   scale, plus vector-difference magnitude and arrows. Save box, resolution
   and mask. Any bubble-relative variant is labelled separately.
6. **Identical geometry:** call `bubble_plot_geometry` once on the shared
   object for 3D and meridional views, annotated with both method labels,
   exact coordinate/mask equality checks, and actual counts. Overlay distinct
   method markers on the same coordinates or state zero maximum coordinate
   difference. Include cap/mobile targets and base/cone sources.

## 8. Validation and implementation stages

1. **Extract shared operators/evaluation.** Record a small legacy matrix,
   right-hand side, solution and validation report before extraction. Compare
   afterward to floating-point tolerance, including existing optional source
   families and legacy scaling settings. No geometry or solver-default change.
2. **Implement explicit assembly.** Assert all block sizes, row slices and
   translation signs, including no cap-normal duplicate and zero translation
   columns in traction/force rows. Test `phi=0`, an interior angle, and `pi`;
   test uniform and cone Stokeslets. Verify `Bm*U` equals `normals(mobile,:)*U`.
3. **Cap-row equivalence.** On deterministic arbitrary `lambda,U`, transform
   Cartesian cap residuals by the exact local `Q_i`; compare each value with
   the legacy split normal/tangential residual before projection. Check equal
   isotropic weighted squared norms and `D_i'*D_i=Pi_i`. Test a weighted case
   that would fail if unequal component weights were mistakenly called equal.
   Also substitute the legacy recovered physical `lambda,U` into the new
   boundary blocks and compare with the legacy boundary residual after its
   documented row permutation/rotation. This checks residual identities,
   not equality of two truncated-SVD minimizers.
4. **Force and field conventions.** Verify `KN'*lambda` sums all base and
   cone coefficients. Integrate a single interior Stokeslet's traction on a
   small accurate sphere quadrature to obtain `-lambda`, checking the sign
   independently of the force map. Compare the shared evaluator with
   [drop_velocity_block](./drop_velocity_block.m) and
   [drop_traction_block](./drop_traction_block.m), including non-unit viscosity,
   shifted centre, multiple batches, and a known one-Stokeslet reference.
   Verify laboratory velocity decays on successively distant sample shells.
5. **Solve and diagnostics.** Compare the explicit scaled solve with a direct
   small pseudoinverse using the same absolute cutoff. Test force-residual
   reporting with an intentionally unbalanced trial vector; no post-solve
   correction. Confirm reported row counts, retained rank and physical
   residuals. Check clean-bubble and fully capped translation against the
   `4*pi*muOut*a` and `6*pi*muOut*a` drag limits on small deterministic cases.
   Establish numerical test bounds by refinement; do not demand interior-edge
   residuals at machine precision.
6. **Comparison and graphics.** Put assertions beside the script's calculations
   and plots for common nonzero scales, lab/relative frame differences,
   pointwise cap vector norm, empty masks, field-difference norms, paired color
   limits, region masks and exterior masking. Test reusable evaluation and
   residual plotting through their own small no-argument checks, including
   hidden figures. Every new function has its own end-of-file no-argument
   self-test; run Code Analyzer and these targeted tests after implementation.
   The experiment-specific assertions execute when the script is explicitly
   run; keep the production script out of automatic test dispatch.
7. **Produce the experiment when explicitly run.** Solve the fixed matched
   pair first. Audit on the second larger/closer grid. Add explicitly configured
   SVD sensitivity cases (`tau/10`, `10*tau`) at identical geometry. A separately
   labelled paired resolution check raises counts through existing geometry
   controls while preserving source type and refinement law; each resolution
   still uses exactly matched sources and targets between formulations.

## 9. Quantitative agreement statement

Report agreement as measured, including failures. A singular-value cutoff is
not a bound on velocity, traction, force error, or discarded field modes.

- For translation report `DeltaU`, both theoretical errors, and the ratio
  `DeltaU/(eU_old+eU_new)` where nonzero. The triangle-inequality check is a
  consistency check, not independent evidence of accuracy. State the actual
  speed errors even when the two translations agree closely.
- For cap relative velocity and mobile normal/traction, compare the difference
  in the corresponding boundary quantities with the sum of their measured
  independent residuals. These necessary bounds do not establish whole-field
  accuracy; cap traction and mobile tangential velocity have no zero target.
- For physical velocity and tangential-traction fields, report the effect of
  validation-grid refinement, matched solve-resolution refinement, and cutoff
  changes. At fixed audit points let `E_old,E_new` be each method's change
  under solve-resolution refinement and `S_old,S_new` its maximum change over
  the cutoff sensitivity cases, with all quantities in the same max-field
  norm and common physical scale. Print
  `Delta/(E_old+E_new+S_old+S_new)` for velocity and traction, with explicit
  zero-denominator handling, alongside the refined cross-method difference.
  Values at most one support an empirical statement that the discrepancy is
  within observed discretization/solver sensitivity. Values above one flag
  unexplained disagreement. These are estimators, not rigorous error bounds;
  stable audit maxima and decreasing differences under refinement are required
  before calling the comparison resolved.
- If refinement/sensitivity checks are omitted, or the edge maxima change by
  more than 10%, label the field-agreement conclusion provisional/unresolved.
  Report the two force residuals and force-weight sensitivity if relevant;
  force-row inclusion alone cannot justify a favorable verdict.

The generated results note must give the numbers behind its conclusion,
sample counts and minimum edge distance, rather than an unsupported statement
that the formulations agree within the SVD tolerance.

## 10. Expected experiment outputs

Under [experiments/results/](./experiments/results/), use the new prefix
`oct29_bubble_cap_formulations` so previous outputs remain available:

- `_summary.csv`: one row per method/configuration/resolution/cutoff, with
  physical parameters, translation/theory errors, boundary/force maxima,
  counts, dimensions, rank, relative/absolute cutoffs, all scalings/weights,
  validation minima and refinement flags.
- `_comparison.csv`: common-grid velocity/traction/translation differences,
  cap/mobile breakdowns, sensitivity estimators, ratios and agreement status.
- `_validation.csv`: point coordinates, masks, edge distance, physical fields
  and unfloored pointwise boundary residuals for both methods; identify the
  primary/fine grids. `_profiles.csv` contains shared polar samples and traces.
- `_data.mat`: both solutions, complete options/settings, immutable shared
  geometry, validation/profile/plane grids, evaluated arrays and figure limits,
  singular diagnostics and reference scales. Permit replotting without solves.
- PNG figures for all six visualization groups, including close edge profiles,
  masked meridional fields and geometry views; keep figure handles available.
- `_results.md`: measured accuracy and agreement statement, method/scaling
  definitions, cutoff/rank diagnostics, minimum edge distance, and limits of
  any omitted or unresolved convergence checks. Keep this numerical report
  separate from the mathematical document.

The production experiment remains an explicitly invoked script; its measured
comparison outputs are generated when it runs.
