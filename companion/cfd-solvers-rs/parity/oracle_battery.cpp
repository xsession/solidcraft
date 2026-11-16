// Test oracle only. Exercises upstream battery/{lithium_ion,dfn,pack,thermal3d}.
#include "cfd/battery/dfn.hpp"
#include "cfd/battery/lithium_ion.hpp"
#include "cfd/battery/pack.hpp"
#include "cfd/battery/thermal3d.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
using namespace cfd::battery;
static bool first = true;
template <class V> static void vec(const std::string& k, const V& v) {
    std::cout << (first ? "" : ",\n") << "\"" << k << "\":[";
    first = false;
    std::size_t i = 0;
    for (auto x : v) {
        double d = static_cast<double>(x);
        std::cout << (i++ ? "," : "");
        if (std::isnan(d)) std::cout << "\"nan\""; else if (std::isinf(d)) std::cout << (d > 0 ? "\"inf\"" : "\"-inf\""); else std::cout << std::setprecision(17) << d;
    }
    std::cout << "]";
}
static std::vector<double> bsr(const BatteryStepResult& s) {
    return {s.time_s, s.current_a, s.terminal_voltage_v, s.open_circuit_voltage_v, s.reaction_overpotential_v, s.electrolyte_overpotential_v, s.ohmic_drop_v,
            s.negative_surface_stoichiometry, s.positive_surface_stoichiometry, s.electrolyte_min_concentration_mol_per_m3, s.electrolyte_max_concentration_mol_per_m3, s.temperature_k, s.heat_generation_w};
}
static std::vector<double> dsr(const DfnStepResult& s) {
    return {s.time_s, s.current_a, s.terminal_voltage_v, s.open_circuit_voltage_v, s.reaction_overpotential_v, s.ohmic_drop_v,
            s.negative_surface_stoichiometry, s.positive_surface_stoichiometry, s.electrolyte_min_concentration_mol_per_m3, s.electrolyte_max_concentration_mol_per_m3, s.temperature_k, s.heat_generation_w};
}
int main() {
    std::cout << "{\n";
    SphericalDiffusionParticle p;
    std::vector<double> pv;
    for (int i = 0; i < 5; ++i) { p.step(2e-6, 10.0); pv.push_back(p.average_stoichiometry()); pv.push_back(p.surface_stoichiometry(2e-6)); }
    vec("particle", pv);
    vec("particle_conc", p.concentrations());
    std::vector<double> ocv;
    for (double x = 0.0; x <= 1.0; x += 0.125) { ocv.push_back(graphite_ocv_v(x)); ocv.push_back(nmc_ocv_v(x)); }
    vec("ocv", ocv);
    vec("bv", std::vector<double>{symmetric_butler_volmer_overpotential_v(3.0, 2.0, 298.15), symmetric_butler_volmer_overpotential_v(-0.5, 1.5, 310.0)});
    auto cfg = default_graphite_nmc_config();
    SingleParticleModel spm(cfg);
    std::vector<double> sv;
    for (int i = 0; i < 6; ++i) { auto s = bsr(spm.step(i < 4 ? 2.0 : -1.0, 5.0)); sv.insert(sv.end(), s.begin(), s.end()); }
    auto d = spm.degradation_state();
    sv.insert(sv.end(), {d.sei_thickness_m, d.plated_lithium_mol, d.active_material_fraction, d.crack_damage});
    vec("spm", sv);
    SingleParticleElectrolyteModel spme(cfg);
    std::vector<double> ev;
    for (int i = 0; i < 6; ++i) { auto s = bsr(spme.step(i < 4 ? 3.0 : -1.5, 4.0)); ev.insert(ev.end(), s.begin(), s.end()); }
    vec("spme", ev);
    vec("spme_ce", spme.electrolyte_concentration());
    vec("spme_phie", spme.electrolyte_potential());
    vec("spme_phis", spme.solid_potential());
    SingleParticleElectrolyteModel spme2(cfg);
    std::vector<DriveCyclePoint> cyc{{10.0, 2.0}, {3.5, 0.0}, {6.0, -1.0}};
    auto dc = simulate_drive_cycle(spme2, cyc, 2.0);
    std::vector<double> dv{static_cast<double>(dc.samples.size()), dc.discharged_capacity_ah, dc.electrical_energy_wh};
    auto last = bsr(dc.samples.back()); dv.insert(dv.end(), last.begin(), last.end());
    vec("drive", dv);
    std::vector<double> fr{0.01, 0.1, 1.0, 10.0, 1000.0}, eis;
    for (auto z : battery_impedance_spectrum(fr)) { eis.push_back(z.real()); eis.push_back(z.imag()); }
    vec("eis", eis);
    ThermalSlab1D slab;
    std::vector<double> q(32, 5e5);
    for (int i = 0; i < 4; ++i) slab.step(0.5, q);
    vec("slab", slab.temperature());
    DoyleFullerNewmanModel dfn(default_graphite_nmc_dfn_config());
    std::vector<double> fv{dfn.total_lithium_mol()};
    for (int i = 0; i < 5; ++i) { auto s = dsr(dfn.step(i < 3 ? 3.0 : -2.0, 5.0)); fv.insert(fv.end(), s.begin(), s.end()); }
    fv.push_back(dfn.total_lithium_mol());
    vec("dfn", fv);
    vec("dfn_ce", dfn.electrolyte_concentration());
    vec("dfn_phie", dfn.electrolyte_potential());
    vec("dfn_phis", dfn.solid_potential());
    vec("dfn_j", dfn.local_reaction_current());
    auto s0 = dsr(dfn.state(1.0)); vec("dfn_state", s0);
    // Pack with an imbalanced cell and active balancing.
    std::vector<LithiumIonCellConfig> cells(3, cfg);
    cells[1].negative.particle.initial_stoichiometry = 0.72;
    cells[1].positive.particle.initial_stoichiometry = 0.45;
    PackControl pc; pc.balancing_current_a = 0.2; pc.balancing_minimum_voltage_v = 3.0;
    SeriesBatteryPack pack(cells, pc);
    std::vector<double> pk;
    for (int i = 0; i < 4; ++i) {
        auto r = pack.step(1.5, 5.0);
        pk.insert(pk.end(), {r.time_s, r.terminal_voltage_v, r.terminal_power_w, r.balancing_heat_w, r.net_terminal_charge_ah, r.terminal_energy_wh, r.balancing_energy_wh});
        for (double b : r.balancing_currents_a) pk.push_back(b);
        for (auto& c : r.cells) pk.push_back(c.terminal_voltage_v);
    }
    vec("pack", pk);
    // Thermal grid: anisotropic, mixed boundaries, hotspot source.
    ThermalGridConfig tg; tg.cells = {4, 3, 2}; tg.spacing_m = {0.01, 0.02, 0.015};
    tg.boundaries[0] = {15.0, 290.0}; tg.boundaries[3] = {40.0, 300.0}; tg.boundaries[5] = {5.0, 280.0};
    std::vector<ThermalVoxel> mats(24);
    for (std::size_t i = 0; i < 24; ++i) mats[i] = {1.5e6 + 1e4 * i, {2.0 + 0.1 * i, 0.5, i % 5 == 0 ? 0.0 : 1.2}};
    ThermalGrid3D grid(tg, mats);
    std::vector<double> qq(24, 0.0); qq[5] = 3.0; qq[17] = 1.0;
    std::vector<double> tr;
    for (int i = 0; i < 3; ++i) { auto r = grid.step(2.0, qq); tr.insert(tr.end(), {r.supplied_power_w, r.outward_cooling_power_w, r.stored_energy_change_j, r.energy_balance_error_j, r.maximum_temperature_k, (double)r.iterations}); }
    vec("thermal", tr);
    vec("thermal_T", grid.temperature());
    // Electrothermal pack: two cells, 2x2x2 voxels each in a 4x2x2 grid.
    ThermalGridConfig eg; eg.cells = {4, 2, 2}; eg.spacing_m = {0.01, 0.01, 0.01}; eg.boundaries[2] = {20.0, 295.0};
    std::vector<LithiumIonCellConfig> ec(2, cfg);
    double cap = cfg.thermal.mass_kg * cfg.thermal.heat_capacity_j_per_kg_k;
    ThermalVoxel v; v.volumetric_heat_capacity_j_per_m3_k = cap / 4.0 / 1e-6; v.conductivity_w_per_m_k = {3.0, 3.0, 3.0};
    std::vector<ThermalVoxel> em(16, v);
    std::vector<std::vector<std::size_t>> regions{{0, 4, 8, 12}, {1, 5, 9, 13}};
    ElectrothermalBatteryPack et(ec, PackControl{}, ThermalGrid3D(eg, em), regions);
    std::vector<double> ev2;
    for (int i = 0; i < 3; ++i) {
        auto r = et.step(4.0, 10.0);
        ev2.insert(ev2.end(), {r.electrical.terminal_voltage_v, r.thermal.supplied_power_w, r.thermal.maximum_temperature_k, r.thermal.energy_balance_error_j, r.electrical.cells[0].temperature_k, r.electrical.cells[1].temperature_k});
    }
    vec("electrothermal", ev2);
    vec("electrothermal_T", et.thermal().temperature());
    std::cout << "\n}\n";
}
