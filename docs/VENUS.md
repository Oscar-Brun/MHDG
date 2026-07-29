# Venus — segregated plasma / neutral coupling

Venus solves the neutral species on its own finite-volume mesh and couples it to
the HDG plasma solver through an outer Picard loop, instead of adding neutral
equations to the monolithic HDG system.

## Why a separate solver

The HDG stabilisation parameter `tau` is built for the strongly anisotropic
plasma transport (parallel >> perpendicular). A neutral fluid is isotropic by
nature. Putting both in one monolithic HDG system forces the neutral equations
through a numerical treatment designed against them.

Venus keeps the plasma side untouched — no plasma equation, boundary condition
or stabilisation is modified — and gives the neutral its own mesh, its own
pseudo-time stepping, its own Riemann solver and its own positivity guarantees.

## Equations

Two neutral models are available, selected by `neutral_model_type`.

**`'Diffusion'`** — the legacy diffusive closure, reproduced on the finite-volume
mesh.

**`'Venus'`** — isothermal Euler, 2 moments. Unknowns are the neutral density
`n_n` and the two poloidal components of the neutral flux `Gamma_n`:

```
d_t n_n     + div(Gamma_n) = n^2 <sv>_rec + S_wall - n_n n <sv>_iz - n_n S_pump

d_t Gamma_n + div(Gamma_n (x) Gamma_n / n_n + c_n^2 n_n I)
            = -nu_fric (Gamma_n - n_n u_i) - nu_iz Gamma_n + n^2 <sv>_rec u_i
```

with the isothermal closure and collision frequencies

```
c_n^2    = Mref * Ti                         (neutrals thermalised on the ions)
nu_fric  = n <sv>_cx + n_n <sv>_nn           (charge exchange + neutral-neutral)
nu_iz    = n <sv>_iz
```

In axisymmetric geometry a `+p_n/R` curvature term is added to the radial
momentum component.

## Relation to the diffusive model

The diffusive model is the exact collisional limit of the Euler model. As
`nu -> infinity` the momentum equation reduces algebraically to

```
Gamma_n -> n_n u_i - grad(c_n^2 n_n) / (nu_fric + nu_iz)
```

and the coefficient `c_n^2 / (nu_fric + nu_iz)` is term for term the `D_nn` of
the diffusive solver, evaluated on the same cell-average state. The Euler model
therefore adds exactly two effects the diffusive closure cannot carry: the drift
`n_n u_i` and the inertia of the neutral flux.

## Numerics

- **Flux**: HLL with Davis wave-speed bounds — native to finite volumes and
  isotropic.
- **Time integration**: local pseudo-time marching to steady state, with a
  per-cell `dt = CFL * V / sum(|S| L)`. Fluxes are explicit; the sources are
  implicit and solved analytically cell by cell, with no Jacobian.
- **Positivity**: the continuity update reads
  `n_new = (n + dt/V (flux + rec + wall)) / (1 + dt/V (iz + pump))`.
  The numerator is positive and the denominator is greater than one, so `n_n > 0`
  holds by construction rather than by clipping.
- **Wall**: mirror state, which makes the HLL mass flux vanish exactly.
  Recycling, puff and pump enter as volume sources, so the mass balance matches
  the diffusive model.

## Coupling

Each time step runs an outer Picard loop: solve the neutral model with the
plasma state frozen, publish `n_n` and the neutral velocity at the plasma
quadrature points, then run the plasma Newton-Raphson with the neutral sources
frozen. The neutral state is under-relaxed by `picard_relax` between iterations.

The plasma reads the neutral field through `compute_neutral_sources`, which is
the same routine the monolithic build uses; `U(phys%idx_rhon_eq)` carries the
neutral density whatever its origin. This is what keeps the monolithic path
unchanged.

## Building and running

Set `MODE = VENUS` and `MDL = $(MDL_NGAMMATITE)` in `lib/Make.inc/arch.make`,
then `make` in `lib/`. See `test/param_venus.txt` for a commented example; the
relevant keys are `neutral_model_type`, `neutral_muscl` and `savePicard` in
`&SWITCH_LST`, and `picard_max_iter`, `picard_tol`, `picard_relax`, `picard_eta`
in `&NUMER_LST`.

Setting `neutral_model_type = 'None'` restores the monolithic legacy behaviour.

## Diagnostics

Each saved solution carries the neutral field under `coupled_neutrals/`
(`rhon`, `Gammax`, `Gammay`). The solver reports the neutral inventory, the wall
fluxes, and the plasma/neutral particle balance, which should close to round-off.

## Status

The finite-volume neutral equation reproduces the monolithic result to 1e-4, and
the plasma/neutral particle balance closes to 1e-14. The Euler model is
first-order in space; the order-2 MUSCL reconstruction with a Venkatakrishnan
limiter is specified but not yet enabled as the default.
