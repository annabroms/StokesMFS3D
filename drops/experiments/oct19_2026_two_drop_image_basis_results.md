# Two-drop image-basis diagnostic — 19 October 2026

Reproduce with `oct19_2026_two_drop_image_basis.m`. MATLAB R2025a;
unit radii, centres `[-1.05 0 0;1.05 0 0]`, active amplitudes one,
active directions `[0 0 1]`, interior viscosity two, exterior viscosity one.
The inner/outer proxy radii are `0.63` and `2`. Four nodes per image line;
requested `N=100` gives 98 proxy nodes per side per drop. There are 182 coarse
collocation points and **1,200 exponentially clustered cap points on 40 rings**
per drop (1,382 total). Independent validation requests 2,000 rotated design
points plus **4,001 Fibonacci gap-cap points** per drop. All seven primary
cases use exactly the same collocation and validation grids, and unweighted
rows. SDR and the joint basis are then checked again with **2,400 cap points
on 57 rings**, keeping the validation grid and source discretisation fixed.

These dense-cap results supersede the provisional 96-cap-point run.

`T_n` is perpendicular to the image line (z for this geometry), `T_t` is
parallel to it (x); both reverse on the opposite drop. Each dense case uses
one full-system SVD with relative cutoff `1e-10`. No one-body SVDs enter
those solves. The scalar comparison below is a maximum over both surfaces;
velocity scales `|S1|a` and traction scale `|S1|` are all one here.

## Independent maximum surface residuals

| Basis / outer-line variant | Continuity | Normal velocity | Tangential traction | Retained rank / columns | Scaled algebraic residual |
|---|---:|---:|---:|---:|---:|
| SDR, mirrored | 9.34593e-03 | 6.48213e-03 | 6.94787e-02 | 1320 / 1320 | 1.15150e-02 |
| SDRT_t, mirrored | 8.35578e-03 | 6.62123e-03 | 6.71503e-02 | 1354 / 1368 | 1.06610e-02 |
| SDRT_n, mirrored | 8.10156e-03 | 6.44339e-03 | 7.09354e-02 | 1362 / 1368 | 1.09540e-02 |
| SDRT_nT_t, mirrored | 7.97404e-03 | 6.20682e-03 | 6.24246e-02 | 1386 / 1416 | 1.02760e-02 |
| SDRT_nT_t, inner images only | 9.54150e-03 | 8.32085e-03 | 6.58057e-02 | 1278 / 1296 | 1.22150e-02 |
| Joint basis, outer endpoint 1.60a | 8.32291e-03 | 6.96686e-03 | 5.82802e-02 | 1380 / 1416 | 1.08620e-02 |
| Joint basis, outer endpoint 1.99a | 7.95521e-03 | 6.20212e-03 | 6.25853e-02 | 1386 / 1416 | 1.02660e-02 |

Unlabelled mirrored endpoints are `1.98a`. All primary matrices have 16,584
rows; the joint mirrored matrix occupies 179.16 MiB, with estimated economy
SVD output storage 209.76 MiB (excluding workspace).

With the dense caps, the joint basis improves the global continuity, normal
velocity, and traction maxima by approximately **15%, 4%, and 10%** relative
to SDR. This revises the provisional traction comparison obtained with only
96 cap points. It remains a modest improvement at this coarse source
resolution. Adding both stresslet directions introduces an exact redundant
mode at every source point; SVD truncation handles it without changing the
physical projection.

## Location of the errors

For every case above, all three global maxima are **outside** the gap-facing
cap of angular radius `pi/5`. Within that cap:

| Basis | Continuity | Normal velocity | Tangential traction |
|---|---:|---:|---:|
| SDR, mirrored | 3.22989e-03 | 3.12413e-03 | 5.29843e-03 |
| Joint, mirrored | 2.60848e-03 | 2.79761e-03 | 4.27730e-03 |

Doubling the collocation cap gives:

| Basis, cap count | Whole continuity | Whole normal | Whole traction | Gap continuity | Gap normal | Gap traction |
|---|---:|---:|---:|---:|---:|---:|
| SDR, 2400 | 9.61839e-03 | 6.75224e-03 | 7.37011e-02 | 2.53620e-03 | 2.53865e-03 | 4.79775e-03 |
| Joint, 2400 | 8.14550e-03 | 6.62326e-03 | 6.33721e-02 | 1.93448e-03 | 2.30854e-03 | 3.88940e-03 |

The joint basis still improves all three global maxima relative to SDR.
However, the residuals are **not invariant under cap refinement**. Unweighted
rows give the denser cap more influence in the least-squares objective, so
this comparison changes both sampling density and the relative weight of
the gap. It cannot by itself establish cap convergence. The largest errors
remain away from the gap even with 2,400 cap points.

This run cannot establish the singularity structure or an accuracy floor at
the demo's higher source resolutions. Coarse proxy and background collocation
refinement is the next useful check, with the dense gap sampling held fixed.
A different fundamental-solution family is not established as necessary by
these results. See the [comparison plot](oct19_2026_two_drop_image_basis_residuals.png).


## Placement and source magnitudes

The first inner radius is `0.7298437881`; its Kelvin image is at
`1.3701562119`. The four outer radii are approximately
`[1.37016, 1.45186, 1.67508, 1.98000]`. The minimum outer clearance from the
owning drop is `0.3701562119`, unchanged by moving the last node. The
geometry audit passed ray alignment, both boundary-side inequalities,
paired counts, Kelvin inversion, monotone radii, and endpoint checks.

Maximum physical coefficient-vector norms for the two mirrored bases:

| Basis / location | Proxy S | Image S | Image D | Image R | Image T_n | Image T_t |
|---|---:|---:|---:|---:|---:|---:|
| SDR inner points | 0.734846 | 1568.22 | 2.51538 | 20.0747 | absent | absent |
| Joint inner points | 0.77968 | 25700.8 | 15.6401 | 115.996 | 13.7001 | 140.951 |
| SDR outer points | 26.1767 | 140.322 | 21.124 | 58.3033 | absent | absent |
| Joint outer points | 44.7236 | 5397.16 | 120.487 | 1425.08 | 167.049 | 484.271 |

The larger coefficients show that a lower least-squares residual can come
with stronger cancellation. Different source families have different
normalisations and dimensions; compare their changes within a family.
Separate directional stresslet coefficients are also nonunique before the
minimum-norm SVD choice. These magnitudes are diagnostics, not errors in
the physical flow.

## Matched iterative check and regression tests

For the default joint mirrored case, one-body GMRES took 19 iterations.
Each local SVD retained 693/708 columns at cutoff `1e-10`. The regularised
GMRES residual was `5.302e-11`; the original scaled algebraic residual was
`1.031e-2`. Independent maxima were `[7.947e-3, 6.126e-3, 6.223e-2]`, and
the absolute Frobenius difference in translations from dense SVD was
`2.671e-5`. Local truncation and the identity-diagonal boundary-data method
need not produce the same least-squares solution as global truncation at
this coarse resolution. A small GMRES residual does not remove
discretisation error.

The deterministic kernel, direction, placement, force-projection, and
matrix-free action tests passed. Integrating each supported kernel's
traction over a sphere verified enclosed-source force/torque normalisations
and zero moments for exterior sources, with maximum error `6.99e-15`.
Both stresslet blocks' shared mode was checked to give zero velocity,
pressure, and traction.
