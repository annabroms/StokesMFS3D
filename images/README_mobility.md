# Image-enhanced mobility for unit spheres

`solve_mobility_with_images(q,Fvec,opt)` computes the six rigid-body
velocities per sphere from prescribed forces and torques. It combines the
recompleted Stokes mobility formulation in [Broms, Barnett and Tornberg,
§§3.2–3.3](https://arxiv.org/pdf/2409.04215) with the close-contact image
discretisation and force extraction in [Broms, Barnett and Tornberg,
§2.3 and Remark 4](https://doi.org/10.1016/j.jcp.2024.113636).
Viscosity and sphere radius are one. Nonunit `opt.viscosity` is rejected by
this mobility solver; the resistance solver and its pressure/traction
evaluators support positive viscosity. `Fvec` has six entries per sphere,
`[Fx;Fy;Fz;Tx;Ty;Tz]`, in the source-strength (force-on-fluid) sign used by
`solve_resistance_with_images`. The returned `U` has
`[vx;vy;vz;omega_x;omega_y;omega_z]` per sphere. Fluid-on-body loads have
the opposite sign.

## Source ordering and moments

The geometry and options are shared with `solve_resistance_with_images`.
Every body starts with coarse proxy Stokeslets. Selected image families
follow in `imageFamilies` order (default `SRD`), with all nodes of one
family contiguous. Each block uses interleaved components
`[x1;y1;z1;x2;y2;z2;...]`. If body `k` has `I_k` image nodes in total and
`b_k` selected families, its vector has `n_k=3(N_proxy+b_k I_k)` entries.
Its `M_k` collocation points give `m_k=3M_k` velocity rows. `useImages=false`
uses only coarse proxy sources and collocation, while still recording close
contacts for validation.

Let `H_k` be the `6-by-n_k` map from raw coefficients to total force and
torque. It includes each basis block's `scale`. For physical Stokeslet force
`f` at `y`, its contribution is `[f;(y-q_k)×f]`; for a rotlet torque `t`, it
is `[0;t]`. Stresslets and potential dipoles have zero columns. The rotlet
kernel is `t×(x-y)/(8*pi*|x-y|^3)` and contributes positive torque in this
convention. This is the image paper's Eq. (28), extended to the ordered
source vector. In particular, the torque constraint combines Stokeslet
moments and rotlet torques.

## Recompletion and rigid motion

Define `C_k=H_k'/(H_k*H_k')` (`n_k-by-6`) and the rank-six projector
`L_k=C_k*H_k` (`n_k-by-n_k`). The implementation applies `L_k*z` as
`C_k*(H_k*z)`, without storing the full projector. With prescribed
`F_k=[force;torque]`, the minimum-norm completion is `lambda0_k=C_k*F_k`.
Thus `H_k*lambda0_k=F_k`, while `H_k*(I-L_k)*z_k=0`. The total physical
strengths and velocity are

```math
\lambda_k=(I-L_k)z_k+\lambda_{0k},\qquad U_k=-H_kz_k.
```

The six components removed from the physical strengths are therefore used
to represent translation and angular velocity. The minus sign and the use
of unprojected `z_k` are mobility-paper Eqs. (49)–(52). Applying `H_k` to
the final physical strengths recovers the prescribed force and torque.
With images disabled, `H_k=getKmat(proxy,q_k)'`, so these are exactly the
basic mobility solver's projection and completion maps.

## Weighted system and matrix-free solve

Let `G_kj` evaluate all sources of body `j` at the collocation points of
body `k`; `K_k=getKmat(collocation_k,q_k)` has size `m_k-by-6`. Let `W_k`
repeat `body.rowScale` for three components at every point. By default,
`opt.useRowScaling=true` sets `body.rowScale=sqrt(body.weights)`;
`false` sets it to one, so `W_k=I`. `rowWeightMode` still selects the
`"paper"` or `"local"` area rule, and the physical `body.weights` remain
stored in the grid in either case. The independent surface check is not
row-scaled. Both mobility demos expose separate resistance and mobility
row-scaling booleans near the top. The
one-body mobility block, of size `m_k-by-n_k`, is

```math
A_k=W_k\{G_{kk}(I-L_k)+K_kH_k\}.
```

`one_body_image_mobility_preconditioner` factors `A_k` with a truncated
SVD. As in the resistance solver, singular values are retained when they
exceed `svdCutoff` times the largest one (default `5e-12`); its factors
apply `A_k^+*mu=Y_k*(U_k'*mu)` without forming a pseudoinverse.
The optionally scaled right-hand side is
`b_k=-W_k*sum_j(G_kj*lambda0_j)`. GMRES applies the square system

```math
\mu_k+\sum_{j\ne k}W_kG_{kj}(I-L_j)A_j^+\mu_j=b_k.
```

The `mu_k` diagonal replaces the rank-deficient rectangular product
`A_k*A_k^+` by identity, following mobility-paper §3.3. Off-diagonal
velocities are summed directly in target batches; the complete global
matrix and pair matrices remain unassembled in the default
`solverMode="oneBody"` path. On convergence, the solver
recovers `z_k=A_k^+*mu_k`, then `U_k` and `lambda_k` by the formulas above.

**Accuracy note.** In the current close-contact examples, this one-body
mobility solve can give a larger residual at independent surface points than
the corresponding image resistance solve. It solves the identity-diagonal
regularized system: a small GMRES residual does not imply an equally small
physical boundary residual. The dense SVD mobility solve below has generally
given higher independent-surface accuracy, though it is much more expensive
in memory and is not guaranteed to improve every discretisation. Compare the
independent `uerr` values when accuracy matters.

For a small dense reference solve, set `solverMode="denseSVD"`. This assembles
the complete row-scaled rectangular system with diagonal blocks `A_k` and
off-diagonal blocks `W_k*G_kj*(I-L_j)`, then applies a truncated global SVD.
Both dense solvers print matrix rows and columns before factorization.
Dense mobility reports one global retained rank, zero GMRES iterations,
`trueGmresResidual=NaN`, and its row-scaled least-squares residual in
`details.linearSystemResidual`.
The independent full-sphere and gap-cap check reports `uerr`; `details`
contains the grid, physical strengths, completion, unprojected strengths,
realised moments, retained ranks, solver residual, and validation samples.
To evaluate flow at new points, use `stokes_mfs_velocity_matrix` with each
body's `details.grid.body(k).basis` and `details.lambda{k}`, or split the
strengths with `image_source_coefficients` and call `direct_image_velocity`.

Compared with `solve_mobility`, this solver handles sphere-specific source
families, contact caps, nonuniform row weights, and unequal one-body blocks.
Compared with `solve_resistance_with_images`, it factors the recompleted
mobility block and adds the six-dimensional completion and velocity map.
Its current path is direct and no-slip; the basic solver remains the path
for ellipsoids, slip, and FMM acceleration.
