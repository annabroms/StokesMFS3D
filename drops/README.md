# Active-drop MFS operators

For near-contact image sources and the new solver interface, see
[README_images.md](README_images.md).

This folder implements a Stokes MFS model for spherical drops with prescribed tangential traction jumps. The exterior fluid has viscosity 1; drop `i` has interior viscosity `mu_i`. Each drop translates but does not rotate. This guide defines the operators as they are implemented here; it is intended to be readable without another document.

## Fields and coefficient ordering

Drop `i` has `N_i` nodes on each proxy surface and `M_i` surface collocation nodes. A vector of nodal 3-vectors is stored with Cartesian components interleaved:

`[v_x(1); v_y(1); v_z(1); v_x(2); v_y(2); v_z(2); ...]`.

The coefficient vector for one drop is `[alpha_i; gamma_i]`, with `3*N_i` entries in each half:

- `alpha_i` are strengths at the inner proxy. They generate the exterior Stokeslet field.
- `gamma_i` are strengths at the outer proxy. They generate the drop's interior field.

For a point `x`, the fields are

- exterior: the sum over drops of the inner-proxy Stokeslet fields with strengths `lambda_i`;
- inside drop `i`: its outer-proxy field with strengths `gamma_i`, divided by `mu_i`, plus the incident exterior fields from every other drop.

The Stokeslet velocity block includes `1/(8*pi)`; `drop_velocity_block` divides by the requested viscosity. The traction block uses the target normal pointing out of the drop. Its kernel has the sign implemented in `drop_traction_block` and does not include a viscosity factor.

## Force-free strength and translation

Let `K_i = kron(ones(N_i,1), eye(3))`, a `3*N_i`-by-3 matrix that repeats a constant vector at every source node. The projected exterior strength is

`lambda_i = Pi_i * alpha_i`, where `Pi_i = I - K_i*K_i'/N_i`.

`Pi_i` removes the three constant-strength modes, so the nodal strengths in `lambda_i` sum to zero. For Stokeslets this makes the exterior representation force-free. It does not impose zero torque.

The removed mean is used for the translation completion:

`U_i = -c_i * K_i' * alpha_i`,

where `c_i` is `completionScale` supplied to the assembly/solver. This couples translation to the otherwise unused constant-strength modes. The corresponding contribution to the normal-velocity residual is `+c_i * normals_i * K_i' * alpha_i`, since the boundary condition is `u_out dot n = U_i dot n`. Spherical experiments commonly set `c_i = 1/(8*pi*a_i)` for radius `a_i`.

## One-drop boundary matrix

For drop `i`, the `6*M_i` rows are grouped as follows:

1. `3*M_i` velocity-continuity components, node-interleaved;
2. `M_i` normal-velocity residuals, one per surface node;
3. `2*M_i` tangential-traction components, with the two tangent directions paired at each node.

The `6*N_i` columns are `[alpha_i; gamma_i]`. Define the following maps, all evaluated at drop `i`'s collocation nodes:

- `G_i`: velocity from the inner proxy with exterior viscosity 1 (`3*M_i` by `3*N_i`);
- `E_i`: velocity from the outer proxy with viscosity `mu_i` (`3*M_i` by `3*N_i`);
- `T_i` and `Q_i`: traction from the inner and outer proxies, respectively (each `3*M_i` by `3*N_i`);
- `B_i`: extracts the normal velocity at each node (`M_i` by `3*M_i`);
- `D_i`: extracts the two tangent components at each node (`2*M_i` by `3*M_i`);
- `C_i = normals_i*K_i'`: maps `alpha_i` to the normal component of its translation-completion term (`M_i` by `3*N_i`).

With these definitions, the self block assembled by `drop_recompleted_matrix` is

$$
A_{ii} =
\begin{bmatrix}
G_i Pi_i & -E_i \\
B_i G_i Pi_i + c_i C_i & 0 \\
D_i T_i Pi_i & -D_i Q_i
\end{bmatrix}.
$$

Thus the first row enforces equality of the two fluid velocities, the second enforces the drop's normal kinematic condition, and the last enforces the tangential traction jump. The prescribed traction data are on the right-hand side, not in this matrix. At node `m`, with active direction `e_i` and amplitude `S_i`, that data is `S_i*(e_i dot n_m)*(e_i dot t_m)` for each of the two tangent directions `t_m`.

The tangent basis is built from the normal in `drop_recompleted_matrix`. The tangential traction is independent of viscosity for a fixed Stokeslet strength; the viscosity contrast appears in the interaction block below.

## Interactions and body ordering

For source drop `j` and target drop `i`, `drop_interaction_block` has the form

$$
A_{ij} =
\begin{bmatrix}
0 & 0 \\
B_i G_{ij} Pi_j & 0 \\
(1-mu_i)D_i T_{ij}Pi_j & 0
\end{bmatrix}, \qquad i \ne j.
$$

Here `G_ij` and `T_ij` evaluate the inner-proxy Stokeslet velocity and traction from drop `j` at drop `i`'s surface. The continuity rows are zero because the incident field from other drops is shared by the inner and outer fluid representations and cancels in the velocity jump. The `gamma_j` columns are zero because only the exterior field of drop `j` acts on another drop. The factor `1-mu_i` is the target drop's viscosity contrast in the traction jump.

`drop_two_body_matrix` assembles `[A_11 A_12; A_21 A_22]`. Rows and columns are body-major: all of body 1's rows/coefficients precede body 2's. It is a dense two-drop reference operator, not the scalable global solver.

## Preconditioned solve

`oneBodyPrecondDrop` row-scales each dense self block before computing its economy SVD. Velocity rows are divided by `velocityScale`; traction rows are multiplied by `tractionRowWeight/tractionScale`. Optional column scaling normalizes nonzero columns before factorization. Singular values are retained when `sigma > svdTolerance*sigma_max`. `applyOneBodyPinv` applies these retained factors in a stable order and undoes column scaling to recover physical `[alpha_i; gamma_i]` coefficients.

The global iteration uses boundary-data vectors `eta_i`, one `6*M_i` block per body. For each matvec, `matvecDropMFS`:

1. applies each local truncated-SVD inverse to `eta_i`;
2. projects the recovered `alpha_i` to `lambda_i = Pi_i*alpha_i`;
3. evaluates only the other-drop velocity and traction contributions;
4. returns `eta` plus those off-diagonal contributions.

The resulting square operator has identity diagonal blocks and cross-drop interactions off diagonal. It avoids using the rank-deficient local map `A_ii*A_ii^+` as the diagonal. `solve_drops_preconditioned` solves this system with GMRES, builds the right-hand side from the prescribed active traction, then recovers `alpha`, `gamma`, `lambda`, and translations. `drop_apply_A` separately applies the row-scaled original operator; its local blocks are cached and its cross-drop action is evaluated directly. These two operators have different roles and should not be conflated.

`drop_validate_surface_conditions` evaluates the physical velocity-continuity, normal-velocity, and tangential-traction residuals on supplied validation points. These are independent of the collocation equations and are the preferred physical accuracy checks.

## Code map

- Kernels: `drop_velocity_block`, `drop_traction_block`.
- Dense assembly: `drop_recompleted_matrix`, `drop_interaction_block`, `drop_two_body_matrix`.
- Local and global solve: `oneBodyPrecondDrop`, `applyOneBodyPinv`, `matvecDropMFS`, `drop_apply_A`, `drop_apply_interactions`, `solve_drops_preconditioned`.
- Field evaluation and validation: `drop_plane_velocity`, `drop_exterior_pressure`, `drop_validate_surface_conditions`.
- Streamlines: `get_streamline`, `compare_streamline_paths`.
- `experiments/` contains single-drop checks, two-drop gap and pressure sweeps, a cluster demo, four-drop flow/streamline studies, and a comparison against an analytic single-drop flow. These scripts are experiments and are not run by the solver functions.
