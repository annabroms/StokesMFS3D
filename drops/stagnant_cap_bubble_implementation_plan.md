# Exterior-only stagnant-cap bubble: implementation plan

This document records the original Stokeslet-only implementation for Section 9 of
`/home/anna/Hämtningar/stagnant_cap_mfs-4.pdf` and the Sadhal–Johnson
stagnant-cap bubble limit. It is separate from the mathematical note. The
existing two-phase drop solver and its experiments remain available.
The later mixed-cone-source extension is summarized at the end and in
`README_stagnant_cap.md`; its basis and matrix sizes supersede the
Stokeslet-only formulas below when extra families are selected.

## Files and reuse

| New file | Responsibility | Existing code reused |
| --- | --- | --- |
| `bubble_options.m` | Validate bubble controls and create a geometry-only adapter | `stagnant_cap_options` |
| `bubble_geometry.m` | Keep the uniform interior proxy, optionally add interior cone rings, and expose point/ring/region counts | `stagnant_cap_geometry`, `stagnant_cap_refined_geometry`, `stagnant_cap_design` |
| `bubble_basis.m` | Build ordered exterior blocks on the cone sites for the later mixed-family extension | `image_family_list`, `image_stresslet_directions`, image-source kernels |
| `bubble_system.m` | Assemble the direct exterior rows and force projection | `stokes_mfs_velocity_matrix`, `stokes_mfs_traction_matrix`, `image_mobility_maps` |
| `bubble_solve.m` | Row-scale, solve by truncated SVD, and recover physical strengths and translation | `getPseudoFactors`, `stagnant_cap_exact_drag` |
| `bubble_validation_grid.m`, `bubble_validate.m`, `bubble_check_validation.m` | Independently sample and check all imposed conditions | `stagnant_cap_validation_grid`, existing batched kernel calls |
| `bubble_plot_geometry.m`, `bubble_plot_residuals.m`, `bubble_plot_strengths.m` | Geometry, boundary, and physical-strength figures | Established edge, marker, and color conventions |

Every new function has a no-argument deterministic self-test. The two experiment
scripts are `experiments/oct24_stagnant_cap_bubble_phi_sweep.m` and
`experiments/oct24_stagnant_cap_bubble_residuals.m`. They write their tables and
figures under `experiments/results/` only when explicitly run.

## Discrete system and data flow

`bubble_geometry` returns only source points with `|y-c|<a`; temporary outer
points made by the shared geometry builders never enter the bubble basis. A
uniform proxy has `N_base` points. When enabled, five default cone levels at
`y_j=c+(a-d_0 q^j)e_r(phi,psi)` add exactly
`N_cone=n_level*n_az_source` geometrical points and `3*N_cone` scalar
unknowns. Geometry reports source ring indices, cap/mobile collocation region
and ring indices, actual counts, and closest source-to-interface distance.

With `N=N_base+N_cone`, `M=M_c+M_m`, and node-interleaved Cartesian vectors,
the three exterior boundary blocks have `M`, `2*M_c`, and `2*M_m` rows. The
projected matrix is `3*M` by `3*N`. `G` includes `1/muOut`; Stokeslet
traction `T` uses the outward normal and is viscosity independent at fixed
physical strength. The force map `K` sums **all** Stokeslet strengths;
`C=K'/(K*K')`, `lambda0=C*Fg`, `lambda=(I-C*K)*alpha+lambda0`, and
`U=-completionScale*K*alpha`. This gives `F_H=-K*lambda` and exact algebraic
force balance before SVD truncation. Translation terms appear in the normal
and cap tangential rows. No interior coefficients, interior kernels, angular
velocity unknown, or derivative source columns are assembled.

The solve prints row blocks, matrix dimensions, estimated dense storage, SVD
cutoff, and retained rank. The analytic reference uses
`stagnant_cap_exact_drag(phi,0,muOut,a)` only for speed comparison and
normalization; the numerical solve never sets the two-phase viscosity to zero.

## Experiment stages and expected output

1. **Uniform-proxy angle sweep.** Keep source points, physical parameters,
   and `svdTolerance=1e-12` fixed; use 13 angles from `0` to `pi`. Update the
   cap/mobile partition and edge collocation only. Export a CSV with source,
   collocation, validation, matrix, rank, speed, force, torque, and four
   boundary maxima columns. Plot signed speed against equation (53), relative
   speed error, and the four normalized boundary maxima. At `0` and `pi`,
   omit edge rings and show absent conditions as `NaN`, not zero.
2. **Cone-source residual experiment.** Start at `phi=pi/4`, `n_level=5`,
   `d0=0.35a`, `q=0.7`, and `n_az_source=64`. Add
   `ringLevels` off-edge collocation rings on **each** side. If omitted for
   cone geometry, `ringLevels` defaults to `2*n_level`; expose azimuthal
   count and angular extent separately.
   The ring angles are `phi ± collocAngularExtent*collocClusterRatio^k`
   (`k=0,...,ringLevels-1`); the ratio is set in the experiment script.
   Export 3D and meridional geometry, four residual maps, and 3D and
   meridional strength maps. Use one strength color scale across base and cone
   nodes: `log10(max(norm(lambda_j)/norm(Fg),1e-16))`. The CSV includes the
   actual equation/unknown counts and nearest source and validation distances.

For the uniform-proxy sweep, the `ringLevels` added rings on each side use
`phi ± edgeBandAngle*collocClusterRatio^k` (`k=0,...,ringLevels-1`). The
default ratio `0.5` preserves the original spacing. A smaller ratio moves
the innermost *added ring* closer to the cap edge without changing the base
collocation nodes or any source point. The nearest collocation node may be a
base-design point, so both actual and added-ring distances are reported.

For both experiments, use a shifted validation grid with at least `5*M` points
and independent refinement towards both sides of the edge. Normalize normal
and cap tangential velocity by `U_ref=abs(Fg_z/zeta_b)` and mobile tangential
traction by `muOut*U_ref/a`. Compare maxima on a second, denser grid with a
closer edge sample. Print the four changes and flag any above 10%; exact force
balance is not counted as a boundary accuracy result.

## Verification

Run the no-argument self-tests for all new functions. In particular, verify
five-level geometry and duplicate exclusion; interleaved block dimensions;
force projection including cone nodes; equivalence of the projected and
explicit boundary residual on a small geometry; endpoint drag coefficients
`4*pi*muOut*a` and `6*pi*muOut*a`; and convergence of independent validation
maxima. Do not execute either experiment as part of the self-test suite.

## Mixed cone source extension

`bubble_basis` places selectable S, R, D, T, T_t, T_n, T_n2, and Q_s
blocks, plus an optional radial or cone-normal S_y derivative, only at the
interior cone nodes. The uniform proxy remains S. With `F` active cone
blocks, the exterior matrix has `3*M` rows and
`3*(N_base+F*N_cone)` columns. The force projection includes every S
block; torque from R or S_y remains a diagnostic.

The experiment default is five cone levels and `ringLevels=10` on each side.
For additional blocks, bubble geometry increases collocation azimuths by
at least `ceil(collocationRatio*(F-1)*N_cone/(2*ringLevels))`, while retaining
source sites and ring offsets. The residual experiment compares cone S with
S+T_t+T_n+T_n2 on common validation grids. The 64-angle sweep retains its
uniform baseline and compares those cone cases at five interior angles.
Strength figures separate coincident families and normalize each by its
physical coefficient scale. Both experiments print pre-assembly dense
storage estimates and report actual ranks, boundary maxima, and torque.
