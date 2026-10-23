# Two-drop image-source comparison at gap 0.1

The executable setup is `compare_drop_image_source_families.m`. Two unit drops
have centres at `x=±1.05`, active amplitude `S1=2`, active directions
`e1=e2=(0,0,1)`, and interior viscosity 10. The comparison keeps mirrored
image lines with 30 nodes each, proxy radii `0.63a` and `2a`, the 1,000-point
exponential cap, and relative SVD cutoff `1e-12`. The requested proxy counts
700 and 1,000 resolve to 686 and 969 spherical-design nodes, respectively.

Fields are compared with the existing `SDRT_nT_t` basis at identical geometry
and independent surface and off-surface targets in the exterior and both
interiors. For each interior, one pressure constant is removed using the
off-surface targets; the same constant shifts normal traction. Velocity is
scaled by `max|S1|a=2`; pressure and full traction by `max|S1|=2`. Each number
below is the maximum across all six target sets. The CSV linked below gives
each fluid and target set separately.

| Requested N | Added images | Max Δvelocity | Max Δpressure | Max Δtraction | Max interface residual | Retained ranks | GMRES |
| ---: | --- | ---: | ---: | ---: | ---: | --- | ---: |
| 700 | none | 0 | 0 | 0 | 8.77e-6 | 4272, 4272 / 5016 | 21 |
| 700 | `T_n2` | 9.47e-7 | 2.29e-4 | 2.29e-4 | 8.61e-6 | 4273, 4273 / 5196 | 21 |
| 700 | `Q_s` | 1.82e-6 | 1.67e-3 | 1.66e-3 | 9.56e-6 | 4278, 4278 / 5196 | 21 |
| 700 | both | 1.78e-6 | 1.63e-3 | 1.62e-3 | 9.31e-6 | 4280, 4279 / 5376 | 21 |
| 1,000 | none | 0 | 0 | 0 | 3.66e-6 | 5563, 5563 / 6714 | 20 |
| 1,000 | `T_n2` | 6.80e-7 | 1.83e-4 | 2.20e-4 | 3.43e-6 | 5564, 5564 / 6894 | 20 |
| 1,000 | `Q_s` | 6.17e-7 | 4.59e-4 | 4.22e-4 | 3.77e-6 | 5515, 5515 / 6894 | 20 |
| 1,000 | both | 7.76e-7 | 4.31e-4 | 4.49e-4 | 3.55e-6 | 5516, 5516 / 7074 | 20 |

The stated field-difference target is `1e-3`. Every N=1,000 case meets it.
At N=700, the two cases containing `Q_s` exceed it on interior **surface**
pressure and traction; their off-surface interior and exterior differences
remain below it. The independent interface residuals stay near `1e-5` at
N=700 and `3e-6` at N=1,000. This is a measured basis sensitivity at the
coarser proxy count, not evidence of a failed quadrupole kernel.

The [source-magnitude plot](results/oct20_2026_image_source_magnitudes.png)
shows the maximum physical vector strength for the coarse Stokeslet proxy and
every selected image type on the inner and mirrored outer lines. The exact
per-region differences, pressure-gauge constants, ranks, source maxima, and
solver diagnostics are in the [MAT results](results/oct20_2026_image_source_comparison.mat);
the [CSV table](results/oct20_2026_image_field_comparison.csv) gives the
pointwise field maxima by fluid and target set. Source strengths, especially
stresslet strengths, are basis-dependent because the three-direction
stresslet has exact zero-field coefficient modes; physical fields are the
appropriate comparison.

Small separate checks cover axial, transverse, and distinct oblique active
directions on the two drops. Analytic `Q_s` fields were checked against source
derivatives, the Stokes PDE, dense/direct actions, rotation covariance, and
force/torque integrals. The three stresslet null vectors at one node give zero
velocity, pressure, and traction and leave a rank-six stresslet field space.
