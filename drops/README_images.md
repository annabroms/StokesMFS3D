# Image-enhanced active spherical drops

`solve_drops_with_images(q,a,S1,e,opt)` returns `[U,lambda,gamma,grid,info]`.
Rows of `q` are centres, `a` holds individual radii, `S1` holds active-stress
amplitudes, and rows of `e` are active directions. The exterior viscosity is
one; `opt.interiorViscosity` may be scalar or body-specific. `grid.body(i).basis`
contains the sources **inside** drop i that represent exterior flow.
`grid.body(i).interiorBasis` contains sources **outside** drop i that represent
its interior flow. `lambda{i}` and `gamma{i}` are their physical strengths,
with `[x;y;z]` interleaved at each node. Every basis begins with a coarse
Stokeslet proxy, followed by selected S/R/T/D image blocks in
`opt.imageFamilies` order. `mirror_ext=false` leaves the outer proxy alone.
Stresslet (`T`) image directions default to the centre-to-neighbour contact
ray. Set `opt.use_orth_dir_images=true` to use a unit direction orthogonal
to that ray: cross it with the Cartesian axis having its smallest absolute
component, then normalise. Opposite contact rays receive opposite directions.
The same choice is used on interior and mirrored exterior image lines;
their point locations are unchanged.

To include both stresslet directions simultaneously, set
`opt.imageFamilies="SDRT_nT_t"` (or `["S","D","R","T_n","T_t"]`).
Here `T_n` is **perpendicular to the image line**, using the deterministic
cross-product rule above; `T_t` is **parallel to the image line**. These
explicit labels ignore `use_orth_dir_images`, which continues to control
plain `T`. In three dimensions the perpendicular plane has two directions;
`T_n` selects one of them, not a full transverse tensor basis.
Both blocks have independent three-component strengths at each image node,
on both source lines when mirrored images are enabled. All kernel evaluators,
force maps, grid arrows, and source-magnitude panels preserve these labels.

The kernel is symmetric in its direction and strength:
`T(r;n)*t = T(r;t)*n`. Thus two perpendicular directional blocks have one
exact shared mode per image node: strengths `[t;-n]` produce zero velocity,
pressure and traction. The six columns span at most five independent fields.
The SVD truncates this null mode; no extra physical projection is introduced.
Individual stresslet coefficient magnitudes therefore depend on the chosen
minimum-norm representation, and cannot alone establish solution accuracy.

For a pair with physical gap `delta`, the near-contact and adaptive resolution
rules use `delta/min(a_i,a_j)`. `image_line_points` locates the unequal-radius
bispherical accumulation points; its original half-Chebyshev nodes run from
the accumulation point toward the inner proxy. An enabled outer line starts
at Kelvin radius `a_i^2/r_i*` and has nodes
`r_k=r_end-(r_end-a_i^2/r_i*)*cos(pi*k/(2*(N_image-1)))`,
where `r_end=r_outer*(1-outerImageOffsetFraction)` and
`k=0,...,N_image-1`. The final node reaches `r_end` when `N_image>=2`;
one-node lines contain only the Kelvin point. The drop-only option
`outerImageOffsetFraction` defaults to `0.01`, so the final node sits 1%
inside the outer proxy radius. Reducing it moves that node closer to the
proxy without moving the inner image line. Invalid proxy/image order is
rejected. The `fibonacci` and `exponential`
contact caps reuse the shared geometry routines. With `useImages=false`,
both the image lines and cap nodes are absent.

## Checking outer-source geometry and accuracy

The top of `demo_two_drops_with_images.m` exposes `capMode`, `capAngles`,
`capPointCounts`, `capExponent`, and `capLayers`. Fibonacci mode uses two
overlapping, contact-centred caps: the two angles and counts correspond to
the wide and narrow caps. Exponential mode uses the first angle and the sum
of the counts on one cap. Its ring angles are
`theta_l = beta*expm1(capExponent*l/capLayers)/expm1(capExponent)`.
Increasing the positive exponent concentrates rings near contact; decreasing
it spreads them farther out. The exponent and ring count are unused in
Fibonacci mode. The demo reports the continuity maximum's angular distance
from contact on each node set, and its independent validation cap follows
the larger configured cap angle. A peak away from contact may require more
points in the wide cap, a wider narrow cap, or weaker exponential clustering;
its location alone does not identify the cause of the residual. Compare
solve-grid and independent-node errors with source and solver settings fixed.
With `useRowScaling=true`, the default row-weight rule is `"paper"` for
Fibonacci and `"local"` for exponential caps. Set `base.rowWeightMode="local"`
to compare both cap modes using the same area-weight construction.

After constructing a grid, run `disp(drop_image_check_geometry(grid,opt))`.
The report gives each pair's inner/outer node counts, first and last radii,
minimum clearance from the owning sphere, minimum outer spacing, Kelvin
product error `abs(r_inner_first*r_outer_first/a^2-1)`, and ray error. It
rejects inconsistent ownership, sources on the wrong side, nonmonotone
lines, wrong endpoints, and mismatched paired counts. With unit radii and
gap `0.1`, the first inner radius is approximately `0.72984` and the outer
Kelvin radius is `1.37016`. With outer proxy ratio `2` and offset `0.01`,
the last outer node is `1.98` when there are at least two nodes.
These checks establish admissible placement; they do not establish optimal
placement for the drop boundary conditions.

The outer sources belong only to their owner's interior representation.
Their field is smooth throughout that drop, so integration of its traction
over the closed sphere gives zero net force and torque by momentum balance
and symmetry of stress. Consequently `gamma` is unprojected. Stokeslets and
rotlets themselves can carry force and torque when enclosed: with this
code's outward target normals, an enclosed Stokeslet gives integrated force
`-f` and torque `-(y-q) cross f`; an enclosed rotlet gives torque `-c`.
The T, T_n, T_t and D kernels have zero integrated force and torque.
`stokes_mfs_traction_matrix()` verifies these facts by sphere quadrature,
including sources outside the integration sphere.

MFS image locations approximate singularities of the analytically continued
solution, not physical singularities within either fluid domain. The
[rigid-sphere image study](https://arxiv.org/abs/2407.03343) found SDR efficient
for its no-slip problems; this does not establish the same choice for coupled
drop conditions. Exact spherical-drop singularity solutions depend on
viscosity contrast and interface conditions, as illustrated by
[Sprenger et al. (2020)](https://www.fuw.edu.pl/~mklis/publications/021.pdf).
Using the same Kelvin accumulation construction for this active-drop
formulation is a placement hypothesis, not an exact drop reflection formula.
Gap-focused sources can help local gradients while coarse proxy resolution
or SVD truncation limits accuracy elsewhere.

Run `experiments/oct19_2026_two_drop_image_basis.m` for a controlled diagnostic.
It compares SDR, SDRT_t, SDRT_n and SDRT_nT_t at identical line locations,
then disables outer images or varies only the outer-line endpoint. Geometry,
forcing, viscosities, proxy surfaces, collocation nodes, validation nodes,
image count, and SVD cutoff otherwise stay fixed. Each case uses one global
dense SVD. The optional `runIterative` check solves the same joint-basis case
with one-body GMRES. The experiment reports all three independent residuals
globally, on the gap-facing cap, and away from it; retained ranks; algebraic
residuals; and proxy/image magnitudes separately by side and family. It also
plots residual locations and joint-basis source strengths. Its default
cap has 1,200 collocation points on 40 exponentially graded rings. Validation
uses an independent 4,001-point Fibonacci cap plus a rotated spherical
design requested at 2,000 points. A final SDR/joint comparison increases the
collocation cap to 2,400 points on 57 rings, with the validation grid fixed.
Since rows are unweighted, increasing cap density also changes its relative
least-squares weight; this is a sensitivity check, not a pure quadrature
convergence test. The default `N=100` is a small source discretisation,
not a converged reference. Increase proxy and
collocation resolution together before drawing conclusions about an accuracy
floor. Extending the outer line alone does not change its first point's
distance from the boundary.

Results from the default small run are recorded in
[`experiments/oct19_2026_two_drop_image_basis_results.md`](experiments/oct19_2026_two_drop_image_basis_results.md).

## Boundary matrix and projection

Let `N_i` be the coarse proxy count, `L_i` the total image nodes over all
near contacts on body i, `f` the number of selected image families, and `M_i`
the collocation count. The one-body matrix has `6*M_i` rows and
`3*(N_i+f*L_i)+3*(N_i+mirror_ext*f*L_i)` columns. Global rows and columns
are body-major. One body's rows are `[3M continuity; M normal velocity;
2M paired tangential traction]`; columns are `[exterior; interior]`.

`H_i` is the 3-by-exterior-column map that sums only **scaled Stokeslet**
strengths, including image Stokeslets. R/T/D have zero columns because their
kernels decay faster than the Stokeslet and carry no net force. With
`C_i=H_i'/(H_i*H_i')`, the physical exterior strengths are
`lambda_i=(I-C_i*H_i)*alpha_i`; `U_i=-completionScale_i*H_i*alpha_i`.
There is no torque projection. The default completion scale is
`1/(8*pi*a_i)`.

Write `G` for exterior-source velocity, `E` for outer-source unit-viscosity
velocity, `T` and `Q` for their physical-strength tractions. With `B` the
normal extractor and `D` the paired tangent extractor, the self block is

```text
[ G P                -E/mu_i ]
[ B G P + c_i n H      0     ]
[ D T P              -D Q    ]
```

For source j and target i, the cross block has zero continuity and interior
columns. Its other rows are `B G_ij P_j` and
`(1-mu_i) D T_ij P_j`. The active traction right-hand side is
`S1_i*(e_i·n)*(e_i·t)` in each tangent direction. Row scaling is optional;
it is disabled by default to match the basic drop solver. The iterative
solver factors only self blocks and applies cross fields in target batches.
The dense global SVD is a two-drop reference. Before each SVD, the solver
prints the matrix dimensions, its actual MATLAB storage in MiB, and the
estimated storage for the three economy SVD outputs. This excludes SVD
workspace and other live arrays. Retained ranks and relative cutoffs are
printed and returned. In `solverMode="denseSVD"`, the two local blocks are
assembled without one-body SVDs; the complete system receives one SVD.
Set `opt.plot_grid=1` to display one drop from the first near-contact pair's
collocation nodes, proxy sources, and image lines before any factorisation
or solve. When `T` is selected, arrows at both sets of image points show
their stresslet directions.

The demos call `drop_image_plot_source_strengths` after solving. It colours
every source point by `log10` of the Euclidean norm of its physical
three-component strength. Proxy sources have separate panels. With mirrored
images, each source family has one interior-image panel and one
exterior-image panel, so both lines are visible; otherwise it has one
interior-image panel. A zero strength is displayed at `log10(realmin)`.

## Flow and checks

`drop_image_field(x,n,region,grid,lambda,gamma,opt)` returns velocity,
pressure, and traction. `region=0` evaluates the global exterior field;
`region=i` evaluates drop i's interior field plus incident exterior fields.
The self interior velocity is divided by `mu_i`; its pressure and traction
use unscaled physical strengths. Incident pressure and traction inside i
are multiplied by `mu_i`. The free-space pressure gauge tends to zero at
infinity. `drop_image_validate_surface_conditions` checks continuity,
normal kinematics, and tangential traction on independent rotated designs
and gap caps. Run `demo_two_drops_with_images.m` for the two discretisations;
set its `runDenseSVD` flag to compare the iterative and dense solves. Run
`demo_P_drops_with_images.m` for a configurable
number of drops. The latter sets `runDenseSVD=false` at the top; set it true
for an optional full-system reference when `P=2`. Its iterative run always
uses the one-body preconditioner, so disabling the dense solve still computes
local one-body SVDs.

`experiments/oct18_2026_two_drop_image_count_sweep.m` keeps the two-drop
geometry, solver settings, and independent validation grid fixed while
sweeping `N_image`. It plots the three maximum scaled surface residuals and
the maximum source-strength norms at inner and outer locations. Proxy
strengths and each S/R/T/D image family are reported separately. The sweep
starts at one image node and fixes `capPointCounts` so its collocation grid
is identical for every case; the zero-image solver would remove those caps.
