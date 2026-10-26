# One spherical drop with a stagnant cap

`solve_stagnant_cap` solves the axis-aligned, fixed-sphere problem in the
stagnant-cap working note (`stagnant_cap_mfs-3.pdf`, equations 20–42 and
Section 8). The cap is
`0 <= theta < phi`, measured from `+z`; the exterior fluid is at rest at
infinity. The cap translates with the drop and has zero **relative** tangential
velocity. No rotation, normal-stress equation, or surfactant-transport solve is
included. The interior and exterior viscosities are both configurable.

## Conical source refinement

`drops/experiments/oct23_stagnant_cap_refinement.m` compares one prescribed-
gravity solve using only uniform proxies with another using the same base
proxies plus source rings on the cone through the cap edge. Both solves use
identical physical parameters, base grids, SVD tolerance, row scaling, and
one larger independent validation grid. The script prints and writes source,
collocation, unknown, SVD, translation, theoretical-speed-error, force,
torque, and boundary-residual statistics. It also exports common-grid
meridional velocity and source-coefficient magnitude plots. Source-strength
figures are saved separately for the uniform and refined cases. Within each
figure, base proxies and each added source family have separate panels for
each flow representation, because different families can occupy the same
nodes. Each panel has its own colour range and shows the logarithm of the
norm of the scaled three-component coefficient. Each case also
exports the four-panel independent surface-residual map from
`plot_stagnant_cap_residuals`, with an automatic color scale for each
boundary-condition panel. The summary CSV reports the maximum normalized
residual for velocity continuity and no penetration over the whole sphere,
cap tangential velocity on the cap, and tangential traction continuity on
the mobile region. Velocity residuals are divided by `Ustar`; traction is
divided by `muOut*Ustar/a`.

`stagnant_cap_refined_geometry(opt)` creates the refined grid. Its defaults
are `n_level=5`, `colloc_factor=2`, `n_az_source=64`, `n_az_colloc=64`,
`sourceClusterRatio=0.7`, and starting distances `0.45a` on each side.
Level `k` is a radial distance `d0*sourceClusterRatio^(k-1)` from the
interface at polar angle `phi`. The last of five rings is `0.108045a`
from the surface. `innerSourceStartDistance` and
`outerSourceStartDistance` are independent. Collocation adds
`ceil(colloc_factor*n_level*n_family)` rings **on each of the cap and mobile
regions**, with offsets starting at `collocAngularExtent=0.2` radians
and multiplied by `collocClusterRatio=0.8` at each level. The edge itself
is excluded. A ring is a polar or radial level, not a single 3D point.
The stagnant-cap experiment scripts expose this count as `ringLevels`, the
number of collocation rings on **each** side of the edge. The refined
experiments convert it to the geometry builder's `colloc_factor` internally.
Here `n_family` counts the active enhancing source families. With the default
single Stokeslet family and azimuthal counts, there are five source rings and 320
added source points on each side, and ten collocation rings and 640
added surface points in each region. If a ring node coincides with a base
node, the whole ring is rotated in azimuth; all counts are preserved.
The geometry's no-argument self-test displays a 3D view and meridional
section and checks these counts.

Added sources default to Stokeslets. `refinementFamilies` uses the ordered
tokens accepted by `image_family_list`; each family adds three coefficients
per source node and the geometry increases collocation azimuthal count if
needed to retain the configured row-to-column oversampling. Non-Stokeslet
families use the shared image kernels but remain experimental for this
stagnant-cap problem; only Stokeslet refinement has been validated here.
To try the three explicit stresslet directions together, set
`refinementFamilies="ST_tT_nT_n2"`. This replaces the legacy `T` with
`T_t`, `T_n`, and `T_n2`; `"STT_tT_nT_n2"` is rejected because it contains
both representations. Each added family increases the coefficient and
collocation counts, so first try a small `n_level` and `n_az_source`.
The force projection includes every exterior-flow Stokeslet, base and added. The
solver computes the source torque as a symmetry diagnostic; this aligned,
nonrotating problem has no torque equation or angular-velocity unknown.
Stokeslet positions and rotlet coefficients contribute to source torque.
No improvement is assumed from the initial source distance: compare the
shared independent residuals and speed error, then vary the distance and
azimuthal resolution to check sensitivity.

## Directional Stokeslet derivatives

`drops/experiments/oct24_stagnant_cap_derivatives.m` compares enhancing
Stokeslets plus radial source derivatives with the same Stokeslets plus
cone-normal source derivatives at **identical** nodes. The option
`derivative_direction` is `"none"` by default, `"radial"` for the physical-drop
normal (e_r), or `"cone_normal"` for the cone normal (e_\theta). Both source
sides use the same signed direction. The geometry returns
`innerDerivativeDirections`, `outerDerivativeDirections`, and corresponding
per-node `DerivativeDistances` arrays. Its self-test displays representative
arrows in 3D and a meridional section.

Each enhancing source point carries an ordinary Stokeslet and an independent
three-component `S_y` derivative coefficient. The derivative acts on the
source coordinate and is multiplied by that node's interface distance after
differentiation; both coefficients have force units. Derivative columns add no
source points and no net force. They add source torque (d_j q_j\times b_j),
whose negative is the reported hydrodynamic torque. The force projection and
completion still use ordinary Stokeslets for prescribed force.

The comparison uses `n_level=5`, two enhancing families, and 20 refined
collocation rings per cap/mobile region. It prints distinct source points,
family and basis-block counts, scalar unknowns, matrix dimensions, speed,
force, torque, and normalized independent residual maxima separately by
boundary condition and region. One independent validation grid serves both
directions and its minimum cap-edge distance is reported. Each direction
produces a four-panel surface residual map on that grid, with matched colour
limits for each condition across cases. Exact zeros are displayed at the
printed per-condition floor. Velocity maps and maxima use (U_*); traction
uses (\mu_{out}U_*/a). A shared meridional plot shows drop-relative
(u_\theta/U_*) and the mobile edge ratio
(u_\theta/[U_*\sqrt{s/a}]), with the exact zero edge limit marked separately.
Source plots split base Stokeslets, enhancing Stokeslets, and derivatives by
flow side. Colours are `log10(max(norm(b)/norm(Fg),1e-16))`, with matched
limits across the two directions for each panel. Tables and figures are
written under `drops/experiments/results/`.

The 24 October 2026 fixed-grid comparison used 1,484 distinct source points,
3,209 collocation points, 6,372 scalar unknowns, and 18,473 independent
validation points. The closest validation point was (4.89\times10^{-5}a)
from the edge. The results in
`experiments/results/oct24_stagnant_cap_derivatives_summary.csv` are:

| direction | speed error | cap/mobile continuity maxima | cap/mobile no-penetration maxima | cap tangential maximum | mobile traction maximum | normalized torque |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| radial | 1.84% | 0.793 / 0.173 | 0.266 / 0.117 | 0.298 | 10.6 | 0.00344 |
| cone-normal | 1.09% | 0.568 / 0.274 | 0.382 / 0.364 | 0.219 | 13.7 | 0.000874 |

These maxima use the scales stated above. Cone-normal derivatives give a
smaller speed error, cap continuity residual, cap tangential residual, torque,
and derivative coefficient magnitude on this grid. Radial derivatives give
smaller mobile continuity, normal, and traction maxima. Both cases retain
large independent traction residuals and very large coefficients. Their
nonzero torque values quantify the unresolved symmetry error in this aligned
problem; no torque equation was added. This experiment does not establish a
uniformly accurate direction. The residual
maps and near-edge rows in the validation CSV should be inspected before
changing the conical geometry in a separate experiment.

## Usage

From the repository root, run `setpath`, then for example:

```matlab
opt = stagnant_cap_options(struct('muIn',2,'muOut',1,'phi',pi/2, ...
    'sourceRequest',100,'ringLevels',4));
geometry = stagnant_cap_geometry(opt);
Fg = (4*pi*opt.radius^3/3)*opt.densityDifference*opt.gravity;
mobility = solve_stagnant_cap(opt,struct('kind','force','Fg',Fg),geometry);
resistance = solve_stagnant_cap(opt,struct('kind','speed','Uz',-1),geometry);
finestOffset = opt.edgeBandAngle*2^(-max(1,opt.ringLevels-1));
validation = stagnant_cap_validation_grid(opt,geometry.collocationCount, ...
    finestOffset);
report = stagnant_cap_validate(mobility,validation,opt);
disp(report.table)
```

The second argument of `stagnant_cap_validation_grid` is the **largest
collocation count**, not the desired validation count; the function selects a
spherical design with at least five times that count when the catalogue has
one, and otherwise uses a deterministic quasiuniform Fibonacci sphere with
the requested number of points. The finer cap-edge validation rings are
added in either case. Run
`drops/experiments/oct23_stagnant_cap.m` for the configured convergence and
limiting-case experiments. It leaves all residual, geometry, and velocity
figures open in MATLAB and also saves PNG files and CSV tables under
`drops/experiments/results/`. Each new function can be called without
arguments to run its own deterministic self-test.

For every resolution, the experiment opens and saves an additional
`oct23_stagnant_cap_collocation_N*.png` figure. Its two 3D views show the
whole physical sphere from opposite sides, with separate markers for base
cap/mobile nodes and added cap/mobile ring nodes.

The convergence experiment also saves one four-panel surface residual map
for each resolution and for each speed/force solve. Each point is an
independent validation point. Colours show `log10` of the pointwise residual
divided by `Ustar` (velocity conditions) or `muOut*Ustar/a` (traction). Each
condition's panel uses its own automatic colour range, so read its colourbar
when comparing panels or resolutions. Exact zeros are drawn one decade below
the smallest positive residual in that panel. The black circle marks
the cap edge. The panels match the note's equations (30)–(33): velocity
continuity and no penetration on the entire sphere, zero exterior tangential
velocity relative to the translating drop on the cap, and zero tangential
traction jump on the mobile interface. The CSV retains the untruncated
statistics and an additional implied interior cap velocity diagnostic, which
is not a separately imposed condition. The script prints one validation table per resolution,
with a `mode` column identifying the prescribed-speed and prescribed-force
rows.

Run `drops/experiments/oct23_stagnant_cap_angle_sweep.m` to compare the
numerical drag force at prescribed `U_z=-1` and the numerical terminal speed
at prescribed gravitational force against equation (50) over cap angle.
Its four panels show each quantity and its absolute relative difference;
the accompanying CSV records actual source and collocation counts and SVD
ranks. Edit `angleSourceRequest`, `angleCount`, or `angleOptions` in that
script to change its resolution or physical parameters. The current script
requests 600 sources and uses an inner proxy radius of `0.6a`.

Run `drops/experiments/oct23_stagnant_cap_N_sweep.m` for a fixed-angle
source-count sweep. Its `useExtraCollocation` boolean at the top controls
whether cap-edge rings are used. The six-panel figure plots, against actual
`N`, the maximum normalized mobile velocity-continuity, no-penetration,
and tangential-traction residuals; the maximum normalized cap exterior
relative tangential velocity; and drag and speed errors against equation
(50). Each boundary panel shows separate prescribed-speed and
prescribed-force curves. One independent validation grid is shared by all
resolutions and retains fine edge sampling even when extra collocation is
off. The result CSV records the six errors, both modes' four residuals,
SVD sizes, and edge distances; a second CSV keeps the near-edge and
elsewhere validation breakdown. Edit `phi`, `sourceRequests`, `ringLevels`,
and physical parameters at the top before running. The main convergence
and angle-sweep scripts also expose `useExtraCollocation` near their tops.

Run `drops/experiments/oct23_stagnant_cap_ring_levels_sweep.m` to hold
`N` and `phi` fixed while varying the number of cap-side and mobile-side
edge-ring levels. Edit `sourceRequest`, `phi`, and `ringLevels` at the top;
the default levels are `1:6`, and level `0` can be included as a base-only
reference. Its eight-panel figure shows the same four boundary maxima and
two equation (50) errors against ring level, plus the total collocation
count `M` and closest collocation distance to the edge. The script checks
that both proxy sets are identical at every level and uses one fixed
independent validation grid. Its CSVs retain counts, SVD sizes, and the
near-edge/elsewhere validation breakdown.

## Exterior-only stagnant-cap bubble experiments

The bubble solver uses only the exterior flow. Uniform proxy sources are
Stokeslets strictly inside the sphere. Optional interior cone nodes carry
the ordered `refinementFamilies` from `image_family_list`: S, R, D, T,
T_t, T_n, T_n2, or Q_s. `derivative_direction` optionally adds an S_y
block in the radial or cone-normal direction. `bubble_options` keeps
`use_orth_dir_images` available for stresslet orientations. Each family
adds three interleaved Cartesian coefficients per cone node; source-site
counts and coefficient counts are reported separately. The force projection
includes every Stokeslet block. Rotlets and derivatives can contribute
source torque, which is reported but not constrained by this nonrotating
formulation.

`ringLevels` is the actual number of collocation rings **in each region**.
At fixed ring offsets, extra cone families raise the azimuthal population
by at least `ceil(collocationRatio*(F-1)*Ncone/(2*ringLevels))` above the
S-only count, where `F` counts cone basis blocks and
`Ncone=n_level*n_az_source`. The geometry also enforces
`3*M >= collocationRatio*coefficientCount`. The default cone comparison
uses five source levels and 10 rings per region. The residual experiment
compares cone S with cone S+T_t+T_n+T_n2 on the same independent validation
grids and exports case-labelled figures plus a comparison CSV. The 64-angle
bubble sweep retains its uniform baseline, adding the two cone cases at
five selected interior angles; endpoints have no cone cases. Both scripts
print dense SVD storage estimates before assembly and must be run explicitly.

`bubble_plot_strengths` keeps base and cone Stokeslets on one colour scale.
Other families use separate panels, with coefficients normalized by
`||F_g|| a` for rotlets and stresslets or `||F_g|| a^2` for dipoles and
quadrupoles. The distance-scaled S_y coefficient uses `||F_g||`.

## Geometry and rows

`get_sphdesign` supplies one uniform spherical design for each concentric
proxy: radius `innerProxyRatio*a` for exterior-flow sources and radius
`outerProxyRatio*a` for interior-flow sources. The base collocation design is
supplemented by geometrically spaced, azimuthally staggered rings on both
sides of the cap edge when `useExtraCollocation=true`. When it is `false`,
only base spherical-design collocation nodes are used. No collocation node
is placed exactly on the edge in either mode.
These earlier experiments do not locally refine source proxies and use only
Stokeslets. The source-refinement comparison described above is separate.
The actual per-surface source count is always `N` (`geometry.N` in MATLAB);
`sourceRequest` is only the requested minimum before the spherical-design
catalogue rounds upward.

With extra collocation on, `M` changes with cap angle. Away from the
endpoints, `M = baseCount + 2*ringLevels*ringAzimuthCount`: every off-edge
base design point is retained, including those inside the edge band, and
the rings add points on both sides. The cap-edge circle has radius
`a*sin(phi)`, so the number of azimuthal nodes per ring rises towards the
equator. At `phi=0` or `pi` there is no interior cap edge, hence no extra
rings. The angle-sweep CSV records all three components of `M`.
With extra collocation off, `M` is the base design count, except that any
base node exactly on the edge is removed.

Let `N` be the source count on each proxy and `Mc`, `Mm` be the cap/mobile
collocation counts, with `M=Mc+Mm`. Nodal vectors use interleaved
`[x1;y1;z1;x2;y2;z2;...]` ordering. `G,T` have size `3M x 3N` and map
exterior Stokeslet strengths to surface velocity and traction. `E,Q` have size
`3M x 3N` and do the same for the interior. `G` includes `1/muOut` and
`E` includes `1/muIn`. For force-normalized Stokeslets, `T,Q` have no
viscosity factor. Both traction blocks contract stress with the normal from
the drop into the exterior. The pointwise tangential projector is
`Pi = I-n*n'`; the two tangent rows of `D` satisfy `D'*D=Pi` per node.
`Nmap` extracts normal velocity, `KM` repeats a translation at interface
nodes, and `KN` repeats a constant vector at inner-proxy source nodes.

The four boundary blocks are velocity continuity on all nodes, zero normal
velocity relative to translation on all nodes, zero cap tangential velocity
relative to translation on cap nodes, and zero tangential traction jump on
mobile nodes. Their row counts are `3M, M, 2Mc, 2Mm`, respectively, for
`6M` rows in total. `stagnant_cap_system` assembles these blocks for both
force and speed experiments. The cap traction jump is an unknown reaction;
the mobile interface can have nonzero common shear traction.

## Force projection and signs

For the force experiment, `Fg=(4*pi*a^3/3)*(rhoIn-rhoOut)*g`,
`PiF=I-KN*KN'/N`, `lambda0=KN*Fg/N`, and
`lambda=PiF*alpha+lambda0`. The drop velocity is
`U=-beta*KN'*alpha`, with default `beta=1/(6*pi*muOut*a)`. This is an
algebraic completion scale, not the drop mobility. The projected normal and
cap tangential rows contain `+beta*Nmap*KM*KN'` and
`+beta*Dc*KM*KN'`, respectively. Every right-hand-side block contains its
known `lambda0` contribution. The resulting matrix is
`6M x 6N`; the force condition is built into `lambda`.

The exterior fluid force on the drop is `FH=-KN'*lambda`; therefore terminal
balance is `FH+Fg=0`. The source-derived torque is
`-sum(cross(yInner-c,lambdaNode))`. For a prescribed signed axial speed,
the same boundary rows solve directly for `[lambda;gamma]`; the measured
resistance is `-FH(3)/Uz`. This independent resistance experiment checks
the gravitational speed. Rows are scaled by `Ustar` for velocity and by
`muOut*Ustar/a` for traction before `getPseudoFactors` applies a truncated
SVD with relative `svdTolerance`.

Before every SVD, `solve_stagnant_cap` prints the scaled matrix's row and
column counts, its real-double storage in decimal MB (`10^6` bytes), the
storage for economy `U`, `S`, and `V`, and a lower bound for known dense
arrays live during the call. That lower bound includes both the assembled
and scaled boundary matrices plus the economy factors; MATLAB/LAPACK SVD
workspace can increase the actual peak. The convergence and angle-sweep
CSVs record the matrix dimensions, input MB, and known-dense lower bound
for each case. Speed and force solves at the same geometry have equal
dimensions and storage estimates, though their retained ranks can differ.

## Accuracy and edge interpretation

`stagnant_cap_validate` reevaluates velocity and traction in batches at an
independent rotated spherical design and finer edge rings. It reports
pointwise maximum and `sqrt(mean(pointwise_norm.^2))`, divided by the stated
velocity or traction scale, separately for cap/mobile and near-edge/elsewhere
subsets. The near-edge subset is within `a*edgeBandAngle`; every row reports
its count and minimum distance from the cap-edge circle. Because rings are
dense, this RMS describes the specified validation sample, not an area-weighted
surface average.

The analytic mobile-side law is `u_theta ~ C*sqrt(s)` for
`s=a*(theta-phi)>0`; its derivative behaves as `C/(2*sqrt(s))`. The surface
velocity tends to zero at the edge. The finite, smooth Stokeslet proxy field
does not reproduce that singular derivative arbitrarily close to the edge.
The velocity figure marks the exact zero separately from the computed trace,
fits `C*sqrt(s)` only over its printed window, and plots `u_theta/sqrt(s)`.
It does not set computed edge values to zero.

Results files produced before the base-node retention change use the old
replacement grid and do not measure the current discretization. Rerun the
experiments to refresh the counts, residuals, and benchmark errors.
