# Image-enhanced rigid-sphere resistance MFS

This folder implements the resistance formulation in Anna Broms, Alex H.
Barnett, and Anna-Karin Tornberg, *Accurate close interactions of Stokes
spheres using lubrication-adapted image systems*, *Journal of Computational
Physics* **523** (2025), 113636. Prescribed
rigid translations and rotations are inputs. Source strengths, forces and
torques are outputs. Viscosity and sphere radius are one.

## Geometry and default discretisation

`getDesignGrid` supplies the unchanged spherical-design proxy sources and
surface collocation nodes. The default `N=686` requests the paper's
686-node coarse design; the corresponding collocation design is selected
by the existing repository routine. The actual proxy count is recorded in
`grid.N`.
The proxy radius is 0.63. For each pair with surface gap `delta < 0.15`, the
default `N_image` is `ceil(-8.72*log10(delta)-6.15)`; otherwise it is zero.
`geometry/spheres/image_line_points.m` implements Eqs. (24), (26), and (27):
nodes lie between `1.05*Rp` and the interior accumulation point and cluster
towards that point. The source families at each image node default to S/R/D
(Stokeslet, rotlet and potential dipole). `imageFamilies` selects a subset of
S/R/T/D; a selected stresslet uses the known line direction as its normal.
Coarse proxy points always carry S. All coefficients are interleaved by
node: `[x1;y1;z1;x2;y2;z2;...]`. Within each body, the order is coarse S
followed by each selected image family in `imageFamilies` order.

Each contact has two caps, with angles `[pi/5,pi/60]` and counts
`[6,6]*3*N_image`, following Section 5.1. `capMode="fibonacci"` reproduces
the paper's quasi-uniform caps. `capMode="exponential"` reuses the earlier
Laplace experiment's graded-ring construction in
`geometry/spheres/collocation_surface_patch.m`. `capAngles`, `capUpsampling`,
`capPointCounts` (an optional exact two-cap count override), `capLayers`,
`capExponent`, `capMinAzimuth`, `N_image`, `nearThreshold`
and `imageFamilies` are editable through the flat solver options.

## Row scaling and one-body solve

For faithful Fibonacci reproduction, `rowWeightMode="paper"` assigns
`4*pi/Mcoarse` to each coarse point and
`2*pi*(1-cos(beta))/Mcap` to each node of a cap. This reproduces the
paper's three area classes per single contact. It is a weighted objective,
not a nonoverlapping surface quadrature: caps can overlap each other and
the coarse grid. Exponential caps must use `rowWeightMode="local"`, which
calls the existing spherical-triangle lumped-area estimator on each body's
complete deduplicated grid (`geometry/spheres/collocation_patch_area.m`).
These positive areas sum to `4*pi` per unit
sphere. The cap comparison applies `"local"` to both cap modes so that
row-weight policy does not confound the node-placement comparison.

The self matrix is `B_k = W_k*G_kk`, with `W_k` repeating `sqrt(area)` for
all three velocity components. A relative singular-value cutoff of
`5e-12` is used by default (Appendix A). `one_body_image_preconditioner`
stores the factored truncated pseudoinverse; the matvec applies
`Y_k*(U_k'*mu_k)` in that order. The preconditioned operator has identity
self blocks and direct-summed off-diagonal interactions, matching Section
2.2. This default `solverMode="oneBody"` path does not assemble the global
matrix or any pair interaction matrix: `image_resistance_matvec` sums S/R/T/D
velocities from source strengths in target batches and reuses pair distances
for image families at the same nodes. It needs no FMM. For a reference solve,
`solverMode="denseSVD"`
assembles the complete weighted matrix and applies a truncated global SVD
using the same `svdCutoff`. It is substantially more costly in memory and
time at the paper's production resolutions.
The image solver prints `starting GMRES`, one dot per completed iteration,
and `completed GMRES` around the iterative solve.

`solve_resistance_with_images(q,U,opt)` returns
`[Fvec,iters,lambda_norm,err_res,details]`. `U` and `Fvec` have six entries
per sphere, `[translation;rotation]` and `[force;torque]` respectively.
`Fvec` follows `solve_resistance` and Eq. (28): sum of Stokeslet forces,
plus source-position moments and rotlet torques. It is the force/torque
applied to the fluid or required to drive the sphere; fluid-on-sphere
force/torque has the opposite sign. T and D carry no net force or torque.
`details` includes body grids, physical source strengths, retained SVD
ranks, the true linear-system residual and independent surface residuals.
`validationVelocity` gives the computed three-component velocity at each
corresponding `validationPoints` row for plotting and further checks.
`iters` is zero for the noniterative dense SVD. Its `retainedRanks` is one
global rank, and `trueGmresResidual` is `NaN`; `linearSystemResidual` is the
true weighted least-squares residual. In one-body mode, `retainedRanks`
contains one rank per sphere and `linearSystemResidual` is the true
preconditioned GMRES residual.
`err_res` follows the paper's pointwise relative velocity residual where
the prescribed speed is nonzero. At a zero-speed validation point, it uses
the largest prescribed speed as a finite denominator.

## Reproduction and checks

Run `demo/demo_resistance_with_images.m` for the public solver example.
At its gap of 0.05, the adaptive rule gives `N_image=6`, hence only 108
nodes per cap under the paper's count formula. For the exponential display,
that demo explicitly requests `[535,543]` cap nodes: the natural counts of
the earlier Laplace 18-ring rule at its two cap angles. Switching the demo
to Fibonacci restores the paper's 108-node caps. The exact-count exponential
routine now preserves `capMinAzimuth` by reducing the layer count if needed.
`demo/demo_cluster_resistance_with_images.m` shows a seeded ten-sphere
cluster with exponential caps and a switch to the paper's Fibonacci caps.
The dated scripts under `experiments/` reproduce the two-sphere cases,
triangle and tetrahedron, and separately compare Fibonacci and exponential
caps on a fixed two-sphere image source geometry. The cap study runs one
proxy-only baseline with no added collocation caps, then compares the two
cap placements with images. Each experiment is an executable MATLAB script
that begins with `clear; close all;`; set its `runFullSweep` flag near the top
to run the larger parameter sweep. The two-sphere resistance and cap-study
scripts also expose `solverMode="oneBody"` or `"denseSVD"` near the top.
Dedicated small tests are
`test_oct16_image_two_body.m` and `test_oct16_image_multibody.m`. New
numerical functions also run deterministic self-tests when called without
arguments. In particular, `solve_resistance_with_images()` prints the
three-sphere forces and surface residuals and plots computed surface
velocities. MATLAB paths are set by the repository's `setpath.m`.
The two-sphere resistance script keeps plotting in its final `%% Plot results`
section. Run that section alone to redraw a `results` struct left
in the workspace or loaded from a MAT file, without repeating the solves.
