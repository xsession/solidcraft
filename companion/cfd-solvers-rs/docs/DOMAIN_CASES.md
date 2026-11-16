# Domain cases (0.7.0)

Version 0.7.0 adds 27 `cfd-solve` case types for the electrical, mechanics, particle,
kinetic and coupled modules ported after 0.6.0. Every one has an embedded example
(`cfd-solve example <name> --out DIR`) and writes `manifest.json` with the
`domain_cases_v1` capability. Inputs are SI unless a field name says otherwise, and the
JSON objects reject unknown fields. Each case exposes one method with narrow inputs; the
library API is broader.

The reference numbers below come from the embedded examples in a release build.
`tests/domain_cases.rs` checks them against closed forms, conservation laws or textbook
values.

## Electrical and electronics

| Solver | Method | Main outputs | Example result |
|---|---|---|---|
| `spice_circuit` | SPICE netlist: `op`, `ac` (log sweep), `transient` (BE, trapezoidal, BDF2), `dc_sweep` | `operating_point.json`, `ac.csv`, `transient.csv`, `dc_sweep.csv` | RC low-pass gain matches `−10 log(1 + (f/fc)²)` within 0.01 dB from 10 Hz to 1 MHz |
| `tcad_diode1d` | 1-D drift–diffusion PN or PIN diode, Gummel iteration, IV and CV sweeps | `equilibrium.csv`, `iv.csv`, `cv.csv` | Built-in potential 0.595 V (2 V_T ln(N/n_i)); forward ideality 1.0–1.3 |
| `battery_cell` | Graphite/NMC SPM, SPMe or DFN cell, or a series SPMe pack; piecewise-constant current with cut-offs | `cycle.csv`, `cycle_summary.json` | 0.25 A discharge delivers 0.286 Ah to 2.8 V |
| `wire_antenna` | Centre-fed straight dipole (or monopole over PEC) by Harrington pulse-basis MoM | `impedance.csv`, `current_distribution.csv` | 0.48 m dipole resonates at 290.8 MHz (0.465 λ) with 71.9 Ω |
| `maxwell3d` | 3-D Yee FDTD, PEC/PMC walls, Gaussian or soft point-pulse excitation, probes | `energy.csv`, `probes.csv`, `fields.vtu` | Field energy is constant to 2 % after the pulse |
| `edge_cavity3d` | Nédélec tetrahedral FEM of a PEC box driven by an edge current port | `port_impedance.csv`, `electric_field.vtu` | First resonance 3.23 GHz on 144 tetrahedra vs TE101 3.35 GHz (−3.6 %) |
| `pdn_ir_drop` | Resistive power-plane grid with VRM sources and loads | `pdn_voltage.vtu` | 11.7 mV maximum drop for 3.5 A |
| `galvanic_corrosion` | Butler–Volmer mixed potential and polarization curves | `polarization.csv`, `mixed_potential.json` | Zn/O₂ couple at −0.476 V, matching the Tafel closed form |
| `power_network` | Per-unit Newton power flow, PV-curve continuation, quadratic dispatch, multi-machine swing | `power_flow.csv`, `pv_curve.csv`, `dispatch.json`, `swing.csv` | Five-bus PV curve reaches 2.8× base load; swing speeds settle below 1e-5 |

`wire_antenna` uses a moment method written for this adapter. The ported
upstream point-matched Pocklington reference (`rf::thin_wire`) remains in the library,
but its feed reactance for a half-wave dipole is tens of kilo-ohms, so it is not exposed.

## Mechanics

| Solver | Method | Main outputs | Example result |
|---|---|---|---|
| `elasticity2d` | Plane stress/strain linear triangles; edge supports, gravity, uniform heating; element stresses | `elasticity.vtu` (displacement, σ, von Mises) | Steel cantilever tip 4.80e-5 m vs Euler–Bernoulli 5.78e-5 m (linear triangles are stiff) |
| `bar_modes1d` | Fixed–free axial bar eigenmodes | `modes.csv`, `mode_shapes.csv` | 1262 Hz vs `c/4L` = 1261.8 Hz |
| `j2_uniaxial` | J2 plasticity with linear hardening under cyclic uniaxial stress (Newton on lateral strain) | `stress_strain.csv` | Elastic to 250 MPa, hardening loops, lateral stress residual < 1e-7 Pa |
| `phase_field_fracture1d` | AT2 phase-field damage in a notched bar under end displacement | `load_displacement.csv`, `damage_profile.csv` | Peak force, then softening to full damage |
| `cosserat_cantilever` | Position-based Cosserat rod clamped over its first segment | `tip_history.csv`, `rod_shape.csv` | Bounded, damped sag; compliance is not calibrated to EI |
| `dem_granular` | Hertz–Mindlin DEM spheres settling in a walled box | `particles_*.vtu`, `particles.pvd`, `settling.csv` | 125 spheres settle; none escape |

## Particles, kinetics and coupled physics

| Solver | Method | Main outputs | Example result |
|---|---|---|---|
| `pic_two_stream1d` | 1-D electrostatic PIC, two cold beams | `energy.csv`, `phase_space.csv`, `fields.csv` | Growth 0.328 ω_p vs cold-beam theory 0.354 ω_p; energy drift 0.3 % |
| `particle_transport_slab` | Condensed-history transport through stacked slabs with sampled discrete interactions | `depth_dose.csv`, `depth_dose.vtu` | Deposited plus secondary energy accounts for the beam within 0.1 % |
| `pennes_bioheat2d` | Pennes bioheat with a Gaussian SAR hot spot; CEM43 dose | `bioheat.vtu`, `hotspot_history.csv` | 40 W/kg for 30 min: +1.2 K peak |
| `joule_heating2d` | DC conduction, Joule source and transient FEM heat | `joule.vtu`, `temperature_history.csv` | Joule power equals `V²σh/w` to 1e-8 |
| `cahn_hilliard1d` | Periodic Cahn–Hilliard spinodal decomposition | `composition.csv`, `free_energy.csv` | Mass conserved; a rising free energy stops the run as unstable |
| `dg_advection1d` | Nodal Gauss–Lobatto DG linear advection | `solution.csv` | Order 4, 16 elements: max error 1.6e-3 after one period |
| `landau_damping1d` | Upwind Vlasov–Poisson, linear Landau damping | `field_energy.csv` | −0.181 vs −0.151 at k = 0.5 (upwind diffusion adds damping) |
| `lennard_jones_md1d` | 1-D Lennard-Jones velocity Verlet with a thermostat phase | `energies.csv`, `final_state.csv` | NVE energy drift 0.16 % |
| `photonic_waveguide` | Scalar finite-difference waveguide modes | `modes.csv`, `modes.vtu` | Si wire n_eff 2.77 on a 16 × 12 grid (coarse; at most 256 cells) |
| `thmc_column1d` | Coupled pressure, temperature, species and reaction column | `history.csv`, `profiles.csv` | Reaction heats the source region |
| `multiphase_box` | Euler–Euler two-fluid column or compressible VOF in a stirred box (FVM hexahedra) | `multiphase_*.vtu`, `multiphase.pvd`, `phase_history.csv` | Phase volume conserved to 1e-15; bubbles rise from 0.09 m to 0.49 m |
| `spray_injection` | Lagrangian spray with drag, heating and d²-law evaporation in a cross-flow | `parcels_*.vtu`, `parcels.pvd`, `vapor_source.vtu`, `spray_history.csv` | 64 % of the injected liquid evaporates; mass balance closes to 1e-13 |

## Limits

Each adapter checks its own size limits (grid cells, particles, steps, sweep points) and
returns an error rather than panicking. None of these cases uses SolidCraft CAD geometry
directly; geometry is generated from the dimensions in the case.

`power_network` accepts an operating point whose residual is below 1e-6 p.u. The ported
finite-difference Newton can stop just above its own 1e-8 target; the strict flag is
reported as `base_newton_strictly_converged`.
