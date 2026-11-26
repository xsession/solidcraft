#!/usr/bin/env python3
"""Writes the benchmark cases in this folder: published problems with known answers, one or
more for each of the solver companion's domains. Every expected value is worked out here, from
the closed form or the published value named with it, so the derivation can be read and run.

    python3 example_cases/benchmarks/make_benchmarks.py
"""
import cmath
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
C0 = 299_792_458.0
Q = 1.602176634e-19
KB = 1.380649e-23
EPS0 = 8.8541878128e-12
BENCHMARKS = []


def benchmark(name, title, domain, shows, source, case, checks):
    BENCHMARKS.append((name, {"title": title, "domain": domain, "shows": shows, "source": source, "case": case, "checks": checks}))


def metric(name):
    return {"metric": name}


def at(file, column, where, value):
    """`column` of `file`, interpolated where `where` is `value`."""
    return {"file": file, "column": column, "at": {"column": where, "value": value}}


def row(file, column, which):
    """`column` of `file` in its `first`, `last`, `max` or `min` row."""
    return {"file": file, "column": column, "row": which}


def near(what, value, expect, tolerance):
    return {"what": what, "value": value, "expect": expect, "tolerance": tolerance}


# --- Structures: axial modes of a fixed-free bar ------------------------------------------
E, RHO, L = 70e9, 2700.0, 0.5
c = math.sqrt(E / RHO)
benchmark(
    "structures-bar-modes", "Axial modes of a fixed–free aluminium bar", "Structures",
    "A 0.5 m aluminium bar fixed at one end vibrates along its length at f_n = (2n − 1)·c/(4L), "
    "c = √(E/ρ): the first four modes on 200 elements.",
    "Fixed–free rod, longitudinal vibration: S. S. Rao, Mechanical Vibrations, 6th ed. (2017), §8.4.",
    {"solver": "bar_modes1d", "length_m": L, "area_m2": 1e-4, "youngs_modulus_pa": E, "density_kg_m3": RHO, "elements": 200, "modes": 4},
    [near(f"mode {n} (Hz)", metric(f"mode_{n}_Hz"), (2 * n - 1) * c / (4 * L), 0.002) for n in range(1, 5)],
)

# --- Plasticity: J2 with linear isotropic hardening, uniaxial ------------------------------
E, SY, H, STRAIN = 200e9, 250e6, 2e9, 0.004
plastic = (STRAIN - SY / E) * E / (E + H)
benchmark(
    "plasticity-j2-uniaxial", "Uniaxial loading past yield with linear hardening", "Plasticity",
    "Steel stretched to 0.4 % strain: elastic to its 250 MPa yield, then σ = σ_y + H·ε_p with "
    "ε_p = (ε − σ_y/E)·E/(E + H).",
    "Return mapping for J2 plasticity with linear hardening: J. C. Simo and T. J. R. Hughes, "
    "Computational Inelasticity (1998), §2.2 (the one-dimensional case).",
    {"solver": "j2_uniaxial", "youngs_modulus_pa": E, "poisson_ratio": 0.3, "yield_stress_pa": SY, "hardening_modulus_pa": H,
     "strain_path": [STRAIN], "increments_per_segment": 80},
    [near("stress at 0.4 % strain (Pa)", metric("max_axial_stress_Pa"), SY + H * plastic, 0.001),
     near("equivalent plastic strain", metric("final_equivalent_plastic_strain"), plastic, 0.001)],
)

# --- Circuits: RC low-pass, frequency response and step response --------------------------
R, CAP = 1e3, 100e-9
fc = 1 / (2 * math.pi * R * CAP)
gain = lambda f: -10 * math.log10(1 + (f / fc) ** 2)
benchmark(
    "circuits-rc-low-pass", "An RC low-pass filter", "Circuits",
    "1 kΩ and 100 nF: the gain falls as −10·log₁₀(1 + (f/f_c)²) with f_c = 1/(2πRC) = 1.59 kHz, "
    "and a 5 V step charges the capacitor to 5·(1 − e^(−t/RC)).",
    "First-order RC circuits: A. S. Sedra and K. C. Smith, Microelectronic Circuits, 7th ed. (2015), §1.6.",
    {"solver": "spice_circuit", "netlist": "* RC low-pass\nV1 in 0 PULSE(0 5 0 1n 1n 1 2) AC 1\nR1 in out 1k\nC1 out 0 100n\n",
     "analyses": [{"type": "ac", "start_hz": 10, "stop_hz": 1e6, "points": 51},
                  {"type": "transient", "dt_s": 1e-7, "steps": 3000, "method": "trapezoidal"}],
     "probes": ["in", "out"]},
    [near("gain at 1 kHz (dB)", at("ac.csv", "out_dB", "frequency_Hz", 1e3), gain(1e3), 0.005),
     near("gain at 10 kHz (dB)", at("ac.csv", "out_dB", "frequency_Hz", 1e4), gain(1e4), 0.005),
     near("output one time constant after the step (V)", at("transient.csv", "V_out", "time_s", R * CAP), 5 * (1 - math.exp(-1)), 0.005)],
)

# --- Semiconductors: an abrupt PN junction -----------------------------------------------
NA, ND, NI, T = 2e21, 1e21, 1e16, 300.0
VT = KB * T / Q
vbi = VT * math.log(NA * ND / NI**2)
current = lambda v: at("iv.csv", "anode_current_A", "anode_voltage_V", v)
benchmark(
    "semiconductors-pn-junction", "An abrupt silicon PN junction", "Semiconductors",
    f"Acceptors 2e21 and donors 1e21 per m³ in silicon (n_i = 1e16 per m³, 300 K): the built-in "
    f"potential is V_T·ln(N_A·N_D/n_i²) = {vbi:.4f} V. Forward biased, the current grows as "
    "exp(V/(n·V_T)) with an ideality n near 1: from 0.3 to 0.4 V by between e^(0.1/(1.1·V_T)) and "
    "e^(0.1/V_T) = 47.9 times.",
    "Abrupt junction and the ideal diode law: S. M. Sze and K. K. Ng, Physics of Semiconductor "
    "Devices, 3rd ed. (2007), §2.2 and §2.3.",
    {"solver": "tcad_diode1d", "length_m": 2e-6, "nodes": 101, "area_m2": 1e-12, "acceptor_m3": NA, "donor_m3": ND,
     "temperature_k": T, "iv": {"start_v": 0.0, "stop_v": 0.5, "points": 11}},
    [near("built-in potential (V)", metric("built_in_potential_V"), vbi, 1e-4),
     {"what": "current at 0.4 V over the current at 0.3 V", "value": {**current(0.4), "over": current(0.3)},
      "at_least": math.exp(0.1 / (1.1 * VT)), "at_most": math.exp(0.1 / VT)}],
)

# --- Antennas: a thin dipole at resonance -----------------------------------------------
L = 0.5
benchmark(
    "antennas-thin-dipole-resonance", "A thin dipole at its resonance", "Antennas",
    "A 0.5 m dipole of 0.1 mm wire (length 5 000 radii) fed at its centre and swept through "
    "half a wavelength: a thin dipole's reactance vanishes when it is about 0.48 wavelengths "
    "long, where its input resistance is about 70 Ω (73 + j42.5 Ω at exactly half a wavelength "
    "is the induced-EMF estimate for a sinusoidal current, which a real wire does not carry).",
    "Thin-dipole resonance: C. A. Balanis, Antenna Theory, 4th ed. (2016), §4.6 and §8.5; "
    "R. W. P. King, The Theory of Linear Antennas (1956).",
    {"solver": "wire_antenna", "length_m": L, "radius_m": 1e-4, "segments": 201,
     "sweep": {"start_hz": 0.9 * C0 / (2 * L), "stop_hz": 1.05 * C0 / (2 * L), "points": 61}, "reference_impedance_ohm": 70},
    [near("length at resonance, in wavelengths", {"metric": "resonant_frequency_Hz", "divide_by": C0 / L}, 0.48, 0.015),
     near("input resistance at resonance (Ω)", metric("resonant_resistance_ohm"), 70.0, 0.05)],
)

# --- EM: the lowest mode of a rectangular cavity ------------------------------------------
A_, B_, D_ = 0.1, 0.08, 0.05
f101 = C0 / 2 * math.sqrt((1 / A_) ** 2 + (1 / D_) ** 2)
benchmark(
    "em-cavity-te101", "The TE101 resonance of a rectangular cavity", "EM",
    "A 100 × 80 × 50 mm box with perfectly conducting walls, driven by a current along its height, "
    f"resonates in TE101 at (c/2)·√((1/a)² + (1/d)²) = {f101 / 1e9:.3f} GHz. The box is cut into "
    "10 mm cubes, each split into six tetrahedra: long, flat boxes converge far more slowly.",
    "Rectangular cavity resonators: D. M. Pozar, Microwave Engineering, 4th ed. (2012), §6.3.",
    {"solver": "edge_cavity3d", "width_m": A_, "height_m": B_, "depth_m": D_, "nx": 10, "ny": 8, "nz": 5,
     "port_a_m": [0.05, 0.04, 0.02], "port_b_m": [0.05, 0.05, 0.02],
     "sweep": {"start_hz": 3.0e9, "stop_hz": 3.6e9, "points": 25}},
    [near("first resonance (Hz)", metric("first_resonance_frequency_Hz"), f101, 0.01)],
)

# --- Compressible flow: Sod's shock tube -------------------------------------------------
def sod_star(gamma=1.4, left=(1.0, 0.0, 1.0), right=(0.125, 0.0, 0.1)):
    """The exact Riemann problem's star state (pressure, velocity, densities), by Newton on p."""
    (rl, ul, pl), (rr, ur, pr) = left, right
    cl, cr = math.sqrt(gamma * pl / rl), math.sqrt(gamma * pr / rr)

    def f(p, r, pk, ck):
        if p > pk:
            a, b = 2 / ((gamma + 1) * r), (gamma - 1) / (gamma + 1) * pk
            return (p - pk) * math.sqrt(a / (p + b))
        return 2 * ck / (gamma - 1) * ((p / pk) ** ((gamma - 1) / (2 * gamma)) - 1)

    p = 0.5 * (pl + pr)
    for _ in range(100):
        g = f(p, rl, pl, cl) + f(p, rr, pr, cr) + ur - ul
        h = 1e-8 * p
        dg = (f(p + h, rl, pl, cl) + f(p + h, rr, pr, cr) - f(p - h, rl, pl, cl) - f(p - h, rr, pr, cr)) / (2 * h)
        p -= g / dg
    u = 0.5 * (ul + ur) + 0.5 * (f(p, rr, pr, cr) - f(p, rl, pl, cl))
    rho_l = rl * (p / pl) ** (1 / gamma)
    g1 = (gamma - 1) / (gamma + 1)
    rho_r = rr * (p / pr + g1) / (g1 * p / pr + 1)
    shock = ur + cr * math.sqrt((gamma + 1) / (2 * gamma) * p / pr + (gamma - 1) / (2 * gamma))
    return p, u, rho_l, rho_r, shock


P, U, RL, RR, S = sod_star()
T_END = 0.2
contact, shock = 0.5 + U * T_END, 0.5 + S * T_END
benchmark(
    "compressible-sod-shock-tube", "Sod's shock tube", "Compressible flow",
    f"A diaphragm at the middle of a tube bursts (ρ, p = 1, 1 left of it and 0.125, 0.1 right, "
    f"γ = 1.4). At t = 0.2 the exact solution has p* = {P:.5f} and u* = {U:.5f} between the "
    f"rarefaction and the shock, density {RL:.5f} left of the contact and {RR:.5f} right of it.",
    "G. A. Sod, J. Comput. Phys. 27 (1978) 1–31; the exact Riemann solution: E. F. Toro, Riemann "
    "Solvers and Numerical Methods for Fluid Dynamics, 3rd ed. (2009), ch. 4.",
    {"solver": "euler1d", "config": {"cells": 400, "reconstruction": "characteristic_weno5", "time_integrator": "ssprk3"},
     "duration": T_END, "left": {"density": 1.0, "velocity": 0.0, "pressure": 1.0}, "right": {"density": 0.125, "velocity": 0.0, "pressure": 0.1},
     "interface_fraction": 0.5},
    [near("pressure between the rarefaction and the shock, x = 0.6", at("flow.csv", "pressure", "x", 0.6), P, 0.01),
     near("density left of the contact, x = 0.6", at("flow.csv", "density", "x", 0.6), RL, 0.01),
     near("density between the contact and the shock", at("flow.csv", "density", "x", 0.5 * (contact + shock)), RR, 0.01)],
)

# --- Kinetic plasma: linear Landau damping -----------------------------------------------
benchmark(
    "kinetic-landau-damping", "Linear Landau damping at k = 0.5", "Kinetic",
    "A small density wave (k·λ_D = 0.5) in a Maxwellian plasma damps without collisions; the "
    "electric field decays at γ = −0.1533 ω_p. The solver's first-order upwind steps add damping "
    "that falls with the cell size: on 512 cells across the wave it is within 4 %.",
    "Root of the Vlasov–Poisson dispersion relation: J. Canosa, J. Comput. Phys. 13 (1973) 158; "
    "used as the standard test in E. Sonnendrücker, Numerical Methods for the Vlasov–Maxwell "
    "Equations (lecture notes, 2015).",
    {"solver": "landau_damping1d", "spatial_cells": 512, "velocity_cells": 128, "wavenumber": 0.5, "amplitude": 0.01, "dt": 0.0025, "steps": 12000},
    [near("damping rate of the field (1/ω_p)", metric("measured_damping_rate"), -0.1533, 0.05)],
)

# --- Plasma: the two-stream instability ---------------------------------------------------
wb = math.sqrt(0.5)
benchmark(
    "plasma-two-stream", "The cold two-stream instability", "Plasma",
    "Two cold electron beams at ±1 m/s (total ω_p = 1/s) in a box whose first mode has "
    "k·v₀ = (√3/2)·ω_b: the fastest-growing wave, growing at ω_b/2 = 0.354 ω_p.",
    "Two-stream dispersion relation: C. K. Birdsall and A. B. Langdon, Plasma Physics via "
    "Computer Simulation (1991), §5.10.",
    {"solver": "pic_two_stream1d", "length_m": 2 * math.pi / (math.sqrt(3) / 2 * wb), "grid_points": 64, "particles_per_beam": 20000,
     "drift_velocity_m_s": 1.0, "plasma_frequency_rad_s": 1.0, "dt_s": 0.05, "steps": 800, "perturbation": 0.001, "mode": 1},
    [near("growth rate (1/s)", metric("measured_growth_rate_per_s"), wb / 2, 0.1)],
)

# --- Lattice Boltzmann: Taylor–Green vortex decay ----------------------------------------
N, TAU, AMP, STEPS = 64, 0.8, 0.02, 400
nu = (TAU - 0.5) / 3
k = 2 * math.pi / N
benchmark(
    "lbm-taylor-green", "Decay of the Taylor–Green vortex", "Flow",
    f"A periodic array of vortices decays by viscosity alone: its kinetic energy falls as "
    f"exp(−4νk²t), ν = (τ − ½)/3 in lattice units, so to {math.exp(-4 * nu * k * k * STEPS):.4f} of its start after {STEPS} steps.",
    "G. I. Taylor and A. E. Green, Proc. R. Soc. Lond. A 158 (1937) 499; for the lattice "
    "Boltzmann method: T. Krüger et al., The Lattice Boltzmann Method (2017), §5.5.",
    {"solver": "lbm_d2q9", "config": {"nx": N, "ny": N, "tau": TAU}, "steps": STEPS, "initial": {"type": "taylor_green", "amplitude": AMP}},
    [near("energy left, as a share of the start", {"metric": "kinetic_energy_lattice", "divide_by": N * N * AMP**2 / 4},
          math.exp(-4 * nu * k * k * STEPS), 0.02)],
)

# --- Chemistry: Robertson's stiff kinetics -----------------------------------------------
def robertson(t_end=40.0):
    from scipy.integrate import solve_ivp

    def f(_, y):
        a, b, c = y
        return [-0.04 * a + 1e4 * b * c, 0.04 * a - 1e4 * b * c - 3e7 * b * b, 3e7 * b * b]

    return solve_ivp(f, (0, t_end), [1.0, 0.0, 0.0], method="Radau", rtol=1e-12, atol=1e-16).y[:, -1]


Y = robertson()
species = [{"name": n, "molar_mass": 0.01, "charge": 0} for n in "ABC"]


def term(i, c=1.0):
    return {"species": i, "coefficient": c}


def rate(k):
    return {"pre_exponential": k, "temperature_exponent": 0.0, "activation_energy": 0.0}


benchmark(
    "chemistry-robertson", "Robertson's stiff reactions", "Chemistry",
    "A → B (0.04/s), B + B → B + C (3e7), B + C → A + C (1e4): rates eleven orders of magnitude "
    f"apart. At t = 40 s: A = {Y[0]:.7f}, B = {Y[1]:.6e}, C = {Y[2]:.7f}, the published reference "
    "(computed again here with an independent stiff solver to 1e-12).",
    "H. H. Robertson, in Numerical Analysis: An Introduction (1966) 178–182; reference values: "
    "E. Hairer and G. Wanner, Solving Ordinary Differential Equations II (1996), and the IVP "
    "test set report ROBER (archimede.uniba.it/~testset).",
    {"solver": "chemistry_batch", "species": species,
     "reactions": [{"reactants": [term(0)], "products": [term(1)], "forward": rate(0.04)},
                   {"reactants": [term(1, 2.0)], "products": [term(1), term(2)], "forward": rate(3e7)},
                   {"reactants": [term(1), term(2)], "products": [term(0), term(2)], "forward": rate(1e4)}],
     "concentrations": [1.0, 0.0, 0.0], "temperature": 300.0, "duration": 40.0,
     "config": {"initial_step": 1e-6, "minimum_step": 1e-16, "relative_tolerance": 1e-7, "absolute_tolerance": 1e-13,
                "max_steps": 20000, "newton_iterations": 8}},
    [near("A at 40 s", {"json": "concentrations.json", "path": ["final", 0]}, Y[0], 1e-4),
     near("B at 40 s", {"json": "concentrations.json", "path": ["final", 1]}, Y[1], 1e-3),
     near("C at 40 s", {"json": "concentrations.json", "path": ["final", 2]}, Y[2], 1e-4)],
)


# --- Optics: spherical aberration of one refracting surface -------------------------------
def trace(h, vertex=0.01, radius=0.05, n1=1.0, n2=1.5, image=0.16):
    """A ray parallel to the axis at height h, refracted by Snell's law at a sphere centred
    on the axis behind its vertex: the ray's height at the plane `image`."""
    centre = vertex + radius
    z = centre - math.sqrt(radius * radius - h * h)
    # The surface normal at the hit, from the centre out: it faces the incoming ray (towards −z).
    nx, nz = h / radius, (z - centre) / radius
    cos_i = -nz
    eta = n1 / n2
    cos_t = math.sqrt(1 - eta * eta * (1 - cos_i * cos_i))
    # Refracted direction: η·d + (η·cos_i − cos_t)·n, with d = (0, 1) and n facing the ray.
    dx, dz = (eta * cos_i - cos_t) * nx, eta + (eta * cos_i - cos_t) * nz
    return h + dx / dz * (image - z)


HEIGHTS = [-4e-3, -2e-3, 0.0, 2e-3, 4e-3]
rms = math.sqrt(sum(trace(h) ** 2 for h in HEIGHTS) / len(HEIGHTS))


def ray(h):
    return {"origin": {"x": h, "y": 0.0, "z": 0.0}, "direction": {"x": 0.0, "y": 0.0, "z": 1.0}, "wavelength_nm": 550.0}


benchmark(
    "optics-spherical-aberration", "Spherical aberration of a refracting sphere", "Optics",
    "Rays parallel to the axis at 0, ±2 and ±4 mm meet a glass sphere of 50 mm radius (n = 1.5); "
    "at the paraxial focus, n₂R/(n₂ − n₁) = 150 mm behind the vertex, the outer rays have "
    f"already crossed the axis: their rms spot radius is {rms * 1e6:.2f} µm, from Snell's law "
    "traced exactly.",
    "Refraction at a spherical surface and spherical aberration: E. Hecht, Optics, 5th ed. "
    "(2017), §5.2 and §6.3.",
    {"solver": "optics_sequential", "object_space_index": 1.0,
     "surfaces": [{"surface_type": "sphere", "vertex_z": 0.01, "radius": 0.05, "aperture_radius": 0.01, "refractive_index_after": 1.5}],
     "rays": [ray(h) for h in HEIGHTS], "image_z": 0.16},
    [near("rms spot radius at the paraxial focus (m)", metric("rms_spot_radius_m"), rms, 0.001)],
)

# --- Fracture: the strength of an AT2 phase-field bar ------------------------------------
E, GC, ELL = 1.0, 1e-3, 0.02
benchmark(
    "fracture-at2-strength", "The strength of a phase-field (AT2) bar", "Fracture",
    "A bar with no notch, pulled at its end: its damage grows evenly, and the stress peaks at "
    "σ_c = (3√3/16)·√(E·G_c/ℓ) before the bar softens.",
    "Homogeneous response of the AT2 model: K. Pham, H. Amor, J.-J. Marigo and C. Maurini, Int. "
    "J. Damage Mech. 20 (2011) 618–652; E. Tanné et al., J. Mech. Phys. Solids 110 (2018) 80–99.",
    {"solver": "phase_field_fracture1d", "cells": 200, "length_m": 1.0, "youngs_modulus_pa": E, "fracture_energy_j_m2": GC,
     "length_scale_m": ELL, "notch_cell": 100, "notch_damage": 0.0, "max_end_displacement_m": 0.25, "load_steps": 1000},
    [near("peak stress (Pa, on a unit section)", metric("peak_force_N"), 3 * math.sqrt(3) / 16 * math.sqrt(E * GC / ELL), 0.002)],
)

# --- Structures: a cantilever under its own weight ----------------------------------------
E, NU, RHO, G, LEN, H = 200e9, 0.3, 7850.0, 9.81, 1.0, 0.1
bending = 1.5 * RHO * G * LEN**4 / (E * H * H)
shear = RHO * G * LEN**2 / (2 * (5 / 6) * E / (2 * (1 + NU)))
benchmark(
    "structures-cantilever-self-weight", "A steel cantilever under its own weight", "Structures",
    "A 1 m × 0.1 m steel plate (plane stress) clamped along one end sags under gravity: "
    f"q·L⁴/(8EI) = {bending * 1e6:.2f} µm from bending and q·L²/(2κGA) = {shear * 1e6:.2f} µm from "
    "shear at its free end.",
    "Cantilever with a uniform load, with the shear term: S. P. Timoshenko and J. M. Gere, "
    "Mechanics of Materials (1972), §9.11.",
    {"solver": "elasticity2d", "width_m": LEN, "height_m": H, "nx": 160, "ny": 16, "youngs_modulus_pa": E, "poisson_ratio": NU,
     "mode": "plane_stress", "density_kg_m3": RHO, "gravity_m_s2": [0, -G], "supports": [{"edge": "left", "ux_m": 0, "uy_m": 0}]},
    [near("largest displacement, at the free end (m)", metric("max_displacement_m"), bending + shear, 0.03)],
)

# --- RF: a mismatched lossy line -------------------------------------------------------
ZC, ZR, VP, ALPHA, LEN = 75.0, 50.0, 2e8, 0.05, 0.03


def s21(f):
    gl = complex(ALPHA, 2 * math.pi * f / VP) * LEN
    a, b, c_, d = cmath.cosh(gl), ZC * cmath.sinh(gl), cmath.sinh(gl) / ZC, cmath.cosh(gl)
    return 2 / (a + b / ZR + c_ * ZR + d)


benchmark(
    "rf-mismatched-line", "A 75 Ω line between 50 Ω ports", "RF",
    "A 30 mm lossy 75 Ω line measured from 50 Ω ports: its transmission S21 = 2/(A + B/Z₀ + C·Z₀ + D) "
    "from the line's ABCD matrix, with the reflections at both ends.",
    "Transmission lines and ABCD parameters: D. M. Pozar, Microwave Engineering, 4th ed. (2012), §2.3 and §4.4.",
    {"solver": "rf_tem", "line": {"characteristic_impedance": ZC, "phase_velocity": VP, "attenuation_nepers_per_m": ALPHA, "length": LEN},
     "frequencies_hz": [1e8, 5e8, 1e9, 2e9], "reference_impedance": ZR},
    [near(f"S21 {part} at {f / 1e6:.0f} MHz", at("s21.csv", f"S21_{part}", "frequency_Hz", f), getattr(s21(f), part), 1e-6)
     for f in (5e8, 2e9) for part in ("real", "imag")],
)

# --- Spray: the d² law -------------------------------------------------------------------
D0, K, T_INJ, T_END = 5e-5, 1e-8, 0.05, 0.15
a_ = K / D0**2
left = ((1 - a_ * (T_END - T_INJ)) ** 2.5 - (1 - a_ * T_END) ** 2.5) / (2.5 * a_ * T_INJ)
benchmark(
    "spray-d2-law", "Droplets evaporating by the d² law", "Spray",
    "50 µm droplets injected for 50 ms into hot air shrink as d² = d₀² − K·t (K = 1e-8 m²/s, a "
    f"250 ms lifetime). 100 ms after the injection stops, {100 * (1 - left):.2f} % of the liquid "
    "has evaporated: the mass of each parcel, (1 − K·t/d₀²)^(3/2), averaged over the injection.",
    "The d² law of droplet evaporation: G. A. E. Godsave, Proc. Combust. Inst. 4 (1953) 818; "
    "S. R. Turns, An Introduction to Combustion, 3rd ed. (2012), ch. 3 and 10.",
    {"solver": "spray_injection", "size_m": [0.4, 0.1, 0.1], "cells": [20, 5, 5], "carrier_velocity_m_s": [2, 0, 0],
     "carrier_temperature_k": 400, "nozzle_m": [0.02, 0.05, 0.05], "injection_velocity_m_s": [10, 0, 0], "droplet_diameter_m": D0,
     "mass_flow_rate_kg_s": 1e-4, "injection_duration_s": T_INJ, "evaporation_constant_m2_s": K, "dt_s": 1e-4,
     "steps": round(T_END / 1e-4), "frames": 5},
    [near("evaporated share of the injected liquid", metric("evaporated_fraction"), 1 - left, 0.01)],
)

# --- Heat: a wall releasing heat, cooled on both faces ------------------------------------
Q_, K_, LEN, RHO, CP = 1e7, 380.0, 0.04, 8900.0, 385.0
VOL = LEN * 0.03 * 0.01
benchmark(
    "heat-wall-with-generation", "A copper wall releasing heat, held at 300 K on both faces", "Heat",
    "10 MW/m³ released evenly in a 40 mm copper wall whose faces are held at 300 K: in the steady "
    "state the temperature is a parabola, q·L²/(8k) = 5.26 K above the faces in the middle and "
    "q·L²/(12k) = 3.51 K above them on average.",
    "Plane wall with uniform generation: T. L. Bergman, A. S. Lavine, F. P. Incropera and D. P. "
    "DeWitt, Fundamentals of Heat and Mass Transfer, 7th ed. (2011), §3.5.",
    {"solver": "solid_heat3d", "mesh": {"dimensions": [40, 1, 1], "lengths": {"x": LEN, "y": 0.03, "z": 0.01}},
     "config": {"dt": 0.05, "density": RHO, "heat_capacity": CP, "conductivity": K_, "temporal_scheme": "backward_bdf2"},
     "output": {"steps": 400, "output_interval": 400}, "initial": 300.0, "source": Q_,
     "boundaries": [{"patch": "left", "kind": "fixed_value", "value": 300.0}, {"patch": "right", "kind": "fixed_value", "value": 300.0}]},
    [near("mean rise above the faces (K)", {"metric": "sensible_energy_J", "add": -RHO * CP * VOL * 300.0, "divide_by": RHO * CP * VOL},
          Q_ * LEN**2 / (12 * K_), 0.01)],
)


# --- Natural convection: de Vahl Davis's square cavity -----------------------------------
def heated_cavity(rayleigh, cells, dt, steps, h=0.1):
    """A square cavity of air (Pr = 0.71), hot on the left, cold on the right, its top and
    bottom adiabatic, gravity down: the temperature difference sets the Rayleigh number. One
    cell deep with slip front and back, so the flow is two-dimensional. The Nusselt number is
    the heat through the hot wall over k·ΔT times the depth."""
    nu, pr, k, rho, g, beta, t_mean = 1.5e-5, 0.71, 0.026, 1.2, 9.81, 1 / 300, 300.0
    alpha = nu / pr
    dt_k = rayleigh * nu * alpha / (g * beta * h**3)
    depth = h / cells
    walls = [{"patch": p, "kind": "fixed_value"} for p in ("left", "right", "bottom", "top")]
    case = {
        "solver": "collocated_flow3d", "method": "piso",
        "mesh": {"dimensions": [cells, cells, 1], "lengths": {"x": h, "y": h, "z": depth}},
        "config": {"dt": dt, "density": rho, "kinematic_viscosity": nu, "pressure_preconditioner": "jacobi",
                   "pressure_iterations": 4000, "pressure_tolerance": 1e-13, "scaled_momentum": True, "piso_neighbour_update": True},
        "output": {"steps": steps, "output_interval": steps},
        "velocity_boundaries": walls + [{"patch": p, "kind": "slip"} for p in ("front", "back")],
        "pressure_boundaries": [],
        "thermal": {"fluid": {"specific_heat_j_per_kg_k": k / (rho * alpha), "conductivity_w_per_m_k": k},
                    "inflow_temperature_k": t_mean, "solids": [],
                    "boundaries": [{"patch": "left", "condition": {"kind": "temperature", "temperature_k": t_mean + dt_k / 2}},
                                   {"patch": "right", "condition": {"kind": "temperature", "temperature_k": t_mean - dt_k / 2}}],
                    "buoyancy": {"gravity_m_per_s2": [0.0, -g, 0.0], "expansion_per_k": beta, "reference_temperature_k": t_mean}},
        "goals": [{"name": "hot", "kind": "patch_heat_out", "target": "left"}, {"name": "cold", "kind": "patch_heat_out", "target": "right"}],
        "goal_control": {"window_steps": 200, "relative_tolerance": 1e-4, "stop_when_converged": True, "min_steps": 1000},
    }
    return case, k * dt_k * depth


for rayleigh, nusselt, cells, dt, steps in [(1e3, 1.118, 32, 0.05, 6000), (1e4, 2.243, 32, 0.05, 6000), (1e5, 4.519, 64, 0.01, 12000)]:
    case, scale = heated_cavity(rayleigh, cells, dt, steps)
    benchmark(
        f"heat-natural-convection-ra{rayleigh:.0e}".replace("+0", ""), f"Natural convection in a square cavity, Ra = {rayleigh:.0e}".replace("+0", ""),
        "Heat",
        f"A square cavity of air, its left wall hot and its right wall cold, the others adiabatic: at "
        f"Rayleigh number {rayleigh:.0e}".replace("+0", "") + f" the air turns in one cell, and the heat crossing it gives a "
        f"mean Nusselt number of {nusselt} ({cells} × {cells} cells, two-dimensional).",
        "G. de Vahl Davis, Natural convection of air in a square cavity: a bench mark numerical "
        "solution, Int. J. Numer. Methods Fluids 3 (1983) 249–264.",
        case,
        [near("mean Nusselt number on the hot wall", {"metric": "goal_hot", "divide_by": -scale}, nusselt, 0.03),
         near("heat out of the cold wall over heat into the hot one", {"metric": "goal_cold", "over": {"metric": "goal_hot"}}, -1.0, 1e-4)],
    )

# --- Acoustics: a Gaussian pressure pulse spreading in two dimensions ----------------------
def pulse_2d(r, t, c, sigma, amplitude):
    """The 2-D wave from p₀ = A·exp(−r²/2σ²) at rest: by the Hankel transform,
    p(r, t) = A·σ²·∫ exp(−k²σ²/2)·cos(c·k·t)·J₀(k·r)·k dk."""
    from scipy.integrate import quad
    from scipy.special import j0

    value, _ = quad(lambda k: math.exp(-0.5 * (k * sigma) ** 2) * math.cos(c * k * t) * j0(k * r) * k, 0, 12 / sigma, limit=400)
    return amplitude * sigma**2 * value


NX, DX, C_SOUND, SIGMA, DT, STEPS = 128, 1e-4, 1500.0, 3e-4, 5e-9, 400
centre = 0.5 * NX * DX
probes = [(30, 0), (32, 0), (34, 0)]  # cells from the centre's lower-left cell, along x and y
checks = []
for i, j in probes:
    x, y = centre + (i + 0.5) * DX, centre + (j + 0.5) * DX
    r = math.hypot(x - centre, y - centre)
    checks.append(near(f"pressure {r * 1e3:.2f} mm from the centre (Pa)", {"vtu": "acoustic", "field": "pressure_Pa", "nearest": [x, y, 0.0]},
                       pulse_2d(r, DT * STEPS, C_SOUND, SIGMA, 1000.0), 0.02))
benchmark(
    "acoustics-gaussian-pulse", "A Gaussian pressure pulse in water", "Acoustics",
    "A 1 kPa Gaussian pulse (σ = 0.3 mm) at rest in water spreads as a ring at 1500 m/s; after 2 µs "
    "the pressure at points on and near the ring is the exact solution of the 2-D wave equation, "
    "from its Hankel transform.",
    "The 2-D initial-value problem: A. D. Pierce, Acoustics (2019), §4.2; the same check is used "
    "for k-space methods by B. E. Treeby and B. T. Cox, J. Biomed. Opt. 15 (2010) 021314.",
    {"solver": "acoustic_kspace2d", "config": {"nx": NX, "ny": NX, "dx": DX, "dy": DX, "dt": DT, "pml_cells": 12},
     "sound_speed": C_SOUND, "density": 1000.0, "steps": STEPS, "center_x_fraction": 0.5, "center_y_fraction": 0.5,
     "sigma": SIGMA, "amplitude": 1000.0, "sensor_cells": [0]},
    checks,
)

# --- Atomistics: a cold Lennard-Jones chain ----------------------------------------------
N_PART, SPACING, T_LOW = 100, 1.12, 0.01
lattice = 4 * sum((n * SPACING) ** -12 - (n * SPACING) ** -6 for n in range(1, 10) if n * SPACING <= 2.5)
benchmark(
    "atomistics-cold-chain", "A cold Lennard-Jones chain", "Atomistics",
    f"100 atoms on a ring at 1.12 σ apart, cooled to T = 0.01 ε/k: each has the lattice energy "
    f"4·Σ((na)⁻¹² − (na)⁻⁶) = {lattice:.5f} ε over its neighbours within the 2.5 σ cutoff, plus T/2 "
    "of potential energy by equipartition in a harmonic crystal.",
    "Lennard-Jones lattice sums and equipartition: N. W. Ashcroft and N. D. Mermin, Solid State "
    "Physics (1976), ch. 20 and 22; M. P. Allen and D. J. Tildesley, Computer Simulation of "
    "Liquids, 2nd ed. (2017), §1.4.",
    {"solver": "lennard_jones_md1d", "particles": N_PART, "box_length": N_PART * SPACING, "initial_temperature": T_LOW,
     "target_temperature": T_LOW, "thermostat_steps": 4000, "dt": 0.005, "steps": 10000},
    [near("potential energy per atom less T/2 (ε)", {"metric": "final_potential_energy_per_particle",
                                                     "minus": {"metric": "mean_temperature", "divide_by": 2.0}}, lattice, 0.002)],
)

# --- Bioheat: uniform heating against perfusion ------------------------------------------
SAR, W_B, RHO_T, C_T, RHO_B, C_B = 40.0, 0.005, 1000.0, 3600.0, 1060.0, 3770.0
rise = RHO_T * SAR / (RHO_B * C_B * W_B)
tau = RHO_T * C_T / (RHO_B * C_B * W_B)
benchmark(
    "bioheat-uniform-heating", "Tissue heated evenly against blood perfusion", "Bioheat",
    f"Tissue absorbing 40 W/kg everywhere, cooled only by perfusion (0.005/s): Pennes's equation "
    f"has it rise as ΔT·(1 − e^(−t/τ)) with ΔT = ρ·SAR/(ρ_b·c_b·ω) = {rise:.3f} K and "
    f"τ = ρc/(ρ_b·c_b·ω) = {tau:.1f} s.",
    "H. H. Pennes, J. Appl. Physiol. 1 (1948) 93–122; the uniform solution of his equation.",
    {"solver": "pennes_bioheat2d", "nx": 16, "ny": 16, "width_m": 0.1, "height_m": 0.1, "peak_sar_w_kg": SAR,
     "hotspot_center_m": [0.05, 0.05], "hotspot_radius_m": 1000.0, "dt_s": 1.0, "steps": 1800, "blood_perfusion_per_s": W_B, "periodic": True},
    [near("rise after one time constant (K)", {**at("hotspot_history.csv", "max_temperature_K", "time_s", tau), "add": -310.15},
          rise * (1 - math.exp(-1)), 0.005),
     near("rise after ten time constants (K)", {"metric": "max_temperature_K", "add": -310.15}, rise * (1 - math.exp(-1800 / tau)), 0.002)],
)

# --- Microstructure: the linear growth of a spinodal wave ---------------------------------
CELLS, KAPPA, A0, DT_CH, STEPS_CH = 32, 5e-4, 1e-4, 1e-4, 3300
k2 = (2 * math.pi) ** 2
sigma_ch = k2 * (0.25 - KAPPA * k2)  # f = c²(1 − c)²/4: f''(½) = −¼
benchmark(
    "microstructure-spinodal-growth", "A spinodal wave growing in a Cahn–Hilliard mixture", "Microstructure",
    "A 50/50 mixture inside its spinodal (f = c²(1 − c)²/4, so f''(½) = −¼) with a small cosine "
    "wave of composition: the wave grows as exp(σt) with σ = M·k²·(−f'' − κ·k²), "
    f"{sigma_ch:.3f} per unit time here.",
    "Linear stability of the Cahn–Hilliard equation: J. W. Cahn, Acta Metall. 9 (1961) 795–801.",
    {"solver": "cahn_hilliard1d", "cells": CELLS, "length": 1.0, "mobility": 1.0, "gradient_energy": KAPPA,
     "mean_composition": 0.5, "perturbation": A0, "dt": DT_CH, "steps": STEPS_CH, "snapshots": 1},
    [near("the wave's amplitude at the end", {"metric": "composition_max", "add": -0.5}, A0 * math.exp(sigma_ch * DT_CH * STEPS_CH), 0.02)],
)

# --- Granular: a sphere resting on a plane, by Hertz's law --------------------------------
R_S, RHO_S, K_N = 0.005, 2500.0, 1e5
mass = RHO_S * 4 / 3 * math.pi * R_S**3
sink = (mass * 9.81 / K_N) ** (2 / 3)
benchmark(
    "granular-hertz-rest", "A sphere resting on a plane", "Granular",
    f"A 10 mm glass bead dropped on the floor settles where its Hertzian contact force k·δ^(3/2) "
    f"carries its weight: it sinks δ = (m·g/k)^(2/3) = {sink * 1e6:.2f} µm.",
    "Hertz's contact law: K. L. Johnson, Contact Mechanics (1985), §4.2.",
    {"solver": "dem_granular", "radius_m": R_S, "density_kg_m3": RHO_S, "count": [1, 1, 1], "container_width_m": 0.05,
     "dt_s": 5e-5, "steps": 20000, "frames": 4, "normal_stiffness": K_N},
    [near("how far it sinks, as a share of Hertz's", {"metric": "final_mean_height_m", "add": -R_S, "divide_by": -sink}, 1.0, 0.02)],
)

# --- Porous media: a slug spreading by diffusion -----------------------------------------
D_S, T_D = 1e-3, 1.0
spread = 2 * math.sqrt(D_S * T_D)
slug = lambda x: 0.5 * (math.erf((x - 0.45) / spread) - math.erf((x - 0.55) / spread))
benchmark(
    "porous-diffusing-slug", "A slug of solute spreading along a column", "Porous media",
    "A slug of solute filling a tenth of a column spreads by diffusion alone (D = 1e-3, no "
    "reaction or coupling): after t = 1 its profile is ½[erf((x − a)/2√(Dt)) − erf((x − b)/2√(Dt))].",
    "J. Crank, The Mathematics of Diffusion, 2nd ed. (1975), §2.2, eq. 2.15.",
    {"solver": "thmc_column1d", "cells": 200, "length": 1.0, "pressure_diffusivity": 1e-3, "thermal_diffusivity": 1e-3,
     "species_diffusivity": D_S, "initial": [1.0, 300.0, 0.0], "source_cells": [90, 109], "source_state": [1.0, 300.0, 1.0],
     "dt": 1e-3, "steps": 1000},
    [near(f"solute at x = {x}", at("profiles.csv", "species", "x", x), slug(x), 0.01) for x in (0.5, 0.56, 0.6)],
)

# --- Turbulence: decaying homogeneous turbulence (k–ε) -------------------------------------
K0, E0, C2, T_TURB = 0.02, 0.01, 1.92, 2.0
k_end = K0 * (1 + (C2 - 1) * E0 * T_TURB / K0) ** (-1 / (C2 - 1))
benchmark(
    "turbulence-k-epsilon-decay", "Decaying turbulence in the k–ε model", "Turbulence",
    "Turbulence with no mean shear only decays: dk/dt = −ε and dε/dt = −C_ε2·ε²/k, so "
    f"k = k₀·(1 + (C_ε2 − 1)·ε₀·t/k₀)^(−1/(C_ε2 − 1)), {k_end / K0:.4f} of its start after 2 s.",
    "B. E. Launder and D. B. Spalding, Comput. Methods Appl. Mech. Eng. 3 (1974) 269–289; S. B. "
    "Pope, Turbulent Flows (2000), §10.4.",
    {"solver": "rans_transport3d", "mesh": {"dimensions": [4, 4, 2], "lengths": {"x": 1.0, "y": 1.0, "z": 0.5}},
     "output": {"steps": 2000, "output_interval": 2000}, "velocity": {"x": 0.0, "y": 0.0, "z": 0.0},
     "velocity_gradient": [0.0] * 9,
     "turbulence": {"model": "k_epsilon", "k": K0, "epsilon": E0, "config": {"transport": {"dt": 1e-3, "temporal_scheme": "crank_nicolson"}}}},
    [near("turbulent kinetic energy after 2 s (m²/s²)", {"vtu": "rans", "field": "k_m2_per_s2", "stat": "mean"}, k_end, 0.01)],
)


# --- Flow: the lid-driven cavity at Re = 100 ---------------------------------------------
GHIA = [(0.8516, 0.23151), (0.6172, -0.13641), (0.5000, -0.20581), (0.4531, -0.21090), (0.2813, -0.15662), (0.1719, -0.10150)]
N_CAV = 65
benchmark(
    "flow-lid-driven-cavity", "The lid-driven cavity at Re = 100", "Flow",
    "A square cavity whose lid slides at 1 m/s (ν = 0.01 m²/s, Re = 100): the air turns in one "
    "large vortex, and its horizontal velocity down the cavity's middle matches Ghia, Ghia and "
    f"Shin's table ({N_CAV} × {N_CAV} cells, two-dimensional). The solver's momentum equation has "
    "first-order upwind convection only, which spreads the vortex: its core turns about 5 % slow "
    "on these cells (within 2 % on 129 × 129, a run of over half an hour).",
    "U. Ghia, K. N. Ghia and C. T. Shin, J. Comput. Phys. 48 (1982) 387–411, Table I.",
    {"solver": "collocated_flow3d", "method": "piso",
     "mesh": {"dimensions": [N_CAV, N_CAV, 1], "lengths": {"x": 1.0, "y": 1.0, "z": 1.0 / N_CAV}},
     "config": {"dt": 0.0075, "kinematic_viscosity": 0.01, "pressure_preconditioner": "jacobi", "pressure_iterations": 4000,
                "pressure_tolerance": 1e-12, "piso_neighbour_update": True},
     "output": {"steps": 8000, "output_interval": 8000},
     "velocity_boundaries": [{"patch": "left", "kind": "fixed_value"}, {"patch": "right", "kind": "fixed_value"},
                             {"patch": "bottom", "kind": "fixed_value"}, {"patch": "top", "kind": "fixed_value", "value": {"x": 1.0, "y": 0.0, "z": 0.0}},
                             {"patch": "front", "kind": "slip"}, {"patch": "back", "kind": "slip"}],
     "goals": [{"name": "mean_speed", "kind": "global_mean_speed"}],
     "goal_control": {"window_steps": 200, "relative_tolerance": 1e-5, "stop_when_converged": True, "min_steps": 1000}},
    [near(f"u at y = {y} on the vertical centreline (m/s)", {"vtu": "flow", "field": "velocity_m_per_s", "component": 0, "nearest": [0.5, y, 0.5 / N_CAV]}, u, 0.06)
     for y, u in GHIA],
)


# --- Power networks: a five-bus load flow --------------------------------------------------
def load_flow(buses, lines, slack_v=1.06):
    """Newton–Raphson load flow in polar form, solved here independently: the voltages."""
    import numpy as np

    n = len(buses)
    y = np.zeros((n, n), complex)
    for a, b, r, x in lines:
        g = 1 / complex(r, x)
        y[a, a] += g
        y[b, b] += g
        y[a, b] -= g
        y[b, a] -= g
    v = np.ones(n, complex)
    v[0] = slack_v
    s_spec = np.array([complex(p, q) for p, q in buses])
    for _ in range(50):
        s = v * np.conj(y @ v)
        mismatch = np.concatenate([(s_spec - s).real[1:], (s_spec - s).imag[1:]])
        if np.max(np.abs(mismatch)) < 1e-12:
            break
        eps = 1e-7
        x0 = np.concatenate([np.angle(v[1:]), np.abs(v[1:])])
        jac = np.zeros((2 * (n - 1), 2 * (n - 1)))

        def residual(xx):
            vv = v.copy()
            vv[1:] = xx[n - 1:] * np.exp(1j * xx[: n - 1])
            ss = vv * np.conj(y @ vv)
            return np.concatenate([(s_spec - ss).real[1:], (s_spec - ss).imag[1:]])

        f0 = residual(x0)
        for col in range(2 * (n - 1)):
            xx = x0.copy()
            xx[col] += eps
            jac[:, col] = (residual(xx) - f0) / eps
        x1 = x0 - np.linalg.solve(jac, f0)
        v[1:] = x1[n - 1:] * np.exp(1j * x1[: n - 1])
    return v


LOADS = [(0.0, 0.0), (-0.45, -0.15), (-0.40, -0.05), (-0.40, -0.05), (-0.60, -0.10)]
LINES = [(0, 1, 0.02, 0.06), (0, 2, 0.08, 0.24), (1, 2, 0.06, 0.18), (1, 3, 0.06, 0.18), (1, 4, 0.04, 0.12), (2, 3, 0.01, 0.03), (3, 4, 0.08, 0.24)]
V = load_flow(LOADS, LINES)
names = ["slack", "north", "lake", "main", "elm"]
benchmark(
    "power-five-bus-load-flow", "A five-bus load flow", "Power networks",
    "The five-bus network of Stagg and El-Abiad (its lines; loads only, no line charging), its "
    "slack bus at 1.06 p.u.: the bus voltages from a Newton–Raphson load flow written here "
    "independently.",
    "Network data: G. W. Stagg and A. H. El-Abiad, Computer Methods in Power System Analysis "
    "(1968); the method: J. J. Grainger and W. D. Stevenson, Power System Analysis (1994), ch. 9.",
    {"solver": "power_network", "buses": [{"name": "slack", "kind": "slack", "voltage_pu": 1.06}]
     + [{"name": nm, "active_pu": p, "reactive_pu": q} for nm, (p, q) in zip(names[1:], LOADS[1:])],
     "lines": [{"from": a, "to": b, "resistance_pu": r, "reactance_pu": x} for a, b, r, x in LINES]},
    [near(f"voltage at {names[i]} (p.u.)", at("power_flow.csv", "voltage_pu", "bus", i), abs(V[i]), 1e-5) for i in range(1, 5)]
    + [near(f"angle at {names[i]} (°)", at("power_flow.csv", "angle_deg", "bus", i), math.degrees(cmath.phase(V[i])), 1e-4) for i in (1, 4)],
)

# --- Batteries: a discharge, then rest -----------------------------------------------------
benchmark(
    "battery-discharge-and-rest", "A cell discharged, then left to rest", "Batteries",
    "A cell discharged at 0.1 A for 15 minutes delivers exactly 0.025 Ah; left to rest for two hours, "
    "its voltage climbs back as the concentrations inside it even out, and settles.",
    "Charge conservation and the relaxation of a lithium-ion cell's overpotentials: M. Doyle, T. "
    "F. Fuller and J. Newman, J. Electrochem. Soc. 140 (1993) 1526; S. G. Marquis et al., J. "
    "Electrochem. Soc. 166 (2019) A3693 (the single-particle model with electrolyte).",
    {"solver": "battery_cell", "model": "spme", "protocol": [{"current_a": 0.1, "duration_s": 900.0}, {"current_a": 0.0, "duration_s": 7200.0}],
     "dt_s": 1.0},
    [near("charge delivered (Ah)", metric("net_discharge_Ah"), 0.025, 1e-6),
     {"what": "recovery in the rest's last half hour (V): it has settled",
      "value": {**row("cycle.csv", "terminal_voltage_V", "last"), "minus": at("cycle.csv", "terminal_voltage_V", "time_s", 900.0 + 7200.0 - 1800.0)},
      "at_least": 0.0, "at_most": 2e-3},
     {"what": "recovery over the whole rest (V): the voltage comes back up",
      "value": {**row("cycle.csv", "terminal_voltage_V", "last"), "minus": at("cycle.csv", "terminal_voltage_V", "time_s", 900.0)},
      "at_least": 1e-3}],
)

# --- High order: discontinuous Galerkin advection ----------------------------------------
benchmark(
    "high-order-dg-advection", "A Gaussian carried once round by DG", "High order",
    "A Gaussian (width 0.08) carried once round a periodic unit domain by order-4 nodal "
    "discontinuous Galerkin on 64 elements comes back to where it began: it matches the exact "
    "translated profile to 2e-6 of its height of 1. Its error falls as h^(p+1): 2^4.9 per halving "
    "from 32 to 64 elements. Mass is kept to rounding.",
    "J. S. Hesthaven and T. Warburton, Nodal Discontinuous Galerkin Methods (2008), §4.8 "
    "(the h^(p+1) error of the upwind scheme for linear advection).",
    {"solver": "dg_advection1d", "elements": 64, "order": 4, "profile": "gaussian", "dt": 2.5e-4, "steps": 4000},
    [{"what": "largest difference from the exact profile after one period", "value": metric("max_error"), "at_most": 2e-6},
     {"what": "mass gained or lost, relative", "value": metric("mass_drift"), "at_least": -1e-12, "at_most": 1e-12}],
)


# --- Photonics: a silicon slab waveguide ---------------------------------------------------
def slab_index(thickness, wavelength=1.55e-6, n_core=math.sqrt(12.11), n_clad=math.sqrt(2.085)):
    """The fundamental TE mode of a symmetric slab: tan(κ·d/2) = γ/κ."""
    from scipy.optimize import brentq

    k0 = 2 * math.pi / wavelength

    def f(n):
        kappa, gamma = k0 * math.sqrt(n_core**2 - n * n), k0 * math.sqrt(n * n - n_clad**2)
        return math.tan(kappa * thickness / 2) - gamma / kappa

    low = math.sqrt(max(n_clad**2, n_core**2 - (math.pi / (k0 * thickness)) ** 2)) + 1e-9
    return brentq(f, low, n_core - 1e-9)


H_WG, W_WG, D_WG, LAMBDA = 1.6e-6, 2e-6, 2.5e-7, 1.55e-6
slab = slab_index(D_WG)
benchmark(
    "photonics-slab-waveguide", "A silicon slab between silica, in a box", "Photonics",
    "A 250 nm silicon layer (n² = 12.11) between silica (n² = 2.085) at 1.55 µm, filling the "
    "width of a 2 µm box whose walls hold the field at zero: the mode separates into the slab's "
    f"TE₀ mode (n = {slab:.4f}, from tan(κd/2) = γ/κ) times sin(πx/W), so "
    f"n_eff = √(n_slab² − (λ/2W)²) = {math.sqrt(slab**2 - (LAMBDA / (2 * W_WG)) ** 2):.4f}. 64 cells across the height.",
    "Symmetric slab waveguide: K. Okamoto, Fundamentals of Optical Waveguides, 2nd ed. (2006), §2.1; "
    "D. Marcuse, Theory of Dielectric Optical Waveguides (1974), ch. 1.",
    {"solver": "photonic_waveguide", "nx": 4, "ny": 64, "width_m": W_WG, "height_m": H_WG,
     "core_m": [0.0, W_WG, (H_WG - D_WG) / 2, (H_WG + D_WG) / 2], "core_permittivity": 12.11, "cladding_permittivity": 2.085,
     "wavelength_m": LAMBDA, "modes": 2},
    [near("fundamental effective index", metric("fundamental_effective_index"), math.sqrt(slab**2 - (LAMBDA / (2 * W_WG)) ** 2), 0.005)],
)

for name, data in BENCHMARKS:
    (HERE / f"{name}.json").write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")
print(f"{len(BENCHMARKS)} benchmarks")
