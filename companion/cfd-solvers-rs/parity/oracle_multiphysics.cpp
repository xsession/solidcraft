// Test oracle only. Never linked into the Rust library or executable.
// Exercises the MIT-licensed upstream multiphysics units and prints JSON goldens.
#include "cfd/fem/mesh2d.hpp"
#include "cfd/fvm/poly_mesh.hpp"
#include "cfd/multiphysics/bioheat.hpp"
#include "cfd/multiphysics/block_system.hpp"
#include "cfd/multiphysics/cht.hpp"
#include "cfd/multiphysics/deformation_optics.hpp"
#include "cfd/multiphysics/electro_thermal.hpp"
#include "cfd/multiphysics/electrochemical_heat.hpp"
#include "cfd/multiphysics/em_spice_bridge.hpp"
#include "cfd/multiphysics/partitioned.hpp"
#include "cfd/multiphysics/porous_electrothermal.hpp"
#include "cfd/multiphysics/projection.hpp"
#include "cfd/multiphysics/transfer.hpp"
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace cfd;
static bool first = true;
static void key(const std::string& k) { std::cout << (first ? "" : ",\n") << "\"" << k << "\":"; first = false; }
static void vec(const std::string& k, const std::vector<double>& v) {
    key(k); std::cout << "[";
    for (std::size_t i = 0; i < v.size(); ++i) std::cout << (i ? "," : "") << std::setprecision(17) << v[i];
    std::cout << "]";
}
static void num(const std::string& k, double v) { key(k); std::cout << std::setprecision(17) << v; }

int main() {
    std::cout << "{\n";
    // Conservative 1-D remap.
    multiphysics::CellGrid1D src{{0.0, 0.3, 0.5, 1.0}}, dst{{0.0, 0.25, 0.6, 0.8, 1.0}};
    vec("transfer_1d", multiphysics::conservative_cell_average_transfer(src, std::vector<double>{1.0, 3.0, 2.0}, dst));
    // Face flux to cell rate.
    auto box = fvm::make_sheared_cartesian_hexa_mesh(3, 2, 1, 0.3, 0.2, 0.1, 0.0);
    std::vector<double> flux(box.face_count());
    for (std::size_t i = 0; i < flux.size(); ++i) flux[i] = 0.1 * static_cast<double>(i) - 0.7;
    vec("face_flux_rate", multiphysics::conservative_face_flux_to_cell_rate(box, flux));
    // Partitioned fixed point: x = cos(x) per component and a coupled linear pair.
    std::vector<double> state{0.5, 1.0, -0.2};
    auto update = [](std::span<const double> x, std::span<double> y) {
        y[0] = std::cos(x[0]);
        y[1] = 0.5 * x[2] + 0.3;
        y[2] = -0.4 * x[1] + 0.1 * std::sin(x[0]);
    };
    auto r = multiphysics::solve_partitioned_fixed_point(state, update);
    vec("partitioned_state", state); num("partitioned_iterations", static_cast<double>(r.iterations)); num("partitioned_relaxation", r.relaxation);
    // Block system.
    multiphysics::BlockCoupledLinearSystem block({2, 1});
    block.add(0, 0, 0, 0, 4.0); block.add(0, 0, 0, 1, -1.0); block.add(0, 1, 0, 0, -1.0); block.add(0, 1, 0, 1, 4.0);
    block.add(0, 1, 1, 0, -1.0); block.add(1, 0, 0, 1, -1.0); block.add(1, 0, 1, 0, 3.0);
    block.add_rhs(0, 0, 1.0); block.add_rhs(0, 1, 2.0); block.add_rhs(1, 0, 3.0);
    std::vector<double> bx(3, 0.0); block.solve(bx); vec("block_solution", bx);
    // Interface laws.
    auto cht = multiphysics::conjugate_interface_heat_flux(350.0, 200.0, 0.002, 300.0, 0.6, 0.001);
    vec("cht", {cht.heat_flux, cht.interface_temperature});
    vec("electrochemical_heat", multiphysics::electrochemical_heat_source(std::vector<double>{10.0, -5.0}, std::vector<double>{0.1, 0.2}, std::vector<double>{1.0, 2.0}));
    multiphysics::PorousElectroThermalFlowState pe{2.0, 300.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i) pe = multiphysics::advance_porous_electrothermal_flow(pe, {}, 0.1, -1000.0);
    vec("porous_electrothermal", {pe.species_concentration_mol_m3, pe.temperature_k, pe.pressure_pa, pe.darcy_velocity_m_per_s});
    auto em = multiphysics::iterate_em_spice_port({0.0, 0.0},
        [](std::complex<double> i) { return std::complex<double>(50.0, 10.0) * i; },
        [](std::complex<double> v) { return (std::complex<double>(5.0, 0.0) - v) / 75.0; });
    vec("em_spice", {em.port_voltage.real(), em.port_voltage.imag(), em.port_current.real(), em.port_current.imag(), static_cast<double>(em.iterations)});
    // Optics deformation.
    optics::SequentialOpticalSystem sys(1.0);
    optics::SequentialSurface s1; s1.type = optics::SurfaceType::sphere; s1.vertex_z = 0.0; s1.radius = 0.05; s1.refractive_index_after = 1.5;
    optics::SequentialSurface s2; s2.type = optics::SurfaceType::plane; s2.vertex_z = 0.01; s2.refractive_index_after = 1.0;
    sys.add_surface(s1); sys.add_surface(s2);
    auto moved = multiphysics::deform_optical_system_axially(sys, [](double z) { return 1e-4 + 0.01 * z; });
    vec("optics_vertices", {moved.surfaces()[0].vertex_z, moved.surfaces()[1].vertex_z});
    // Projections.
    multiphysics::StructuredScalarGrid3D grid{3, 3, 2, 1.0, 1.0, 1.0, {}};
    for (std::size_t i = 0; i < 18; ++i) grid.values.push_back(std::sin(0.7 * static_cast<double>(i)));
    vec("structured_samples", {grid.sample(0.3, 0.6, 0.25), grid.sample(1.2, -0.1, 0.9)});
    auto cells = fvm::make_sheared_cartesian_hexa_mesh(4, 3, 1, 1.0, 1.0, 0.1, 0.0);
    vec("structured_to_cells", multiphysics::structured_to_fvm_cells(grid, cells));
    auto tri = fem::make_rectangle_tri_mesh(4, 4, 1.0, 1.0);
    std::vector<double> nodal(tri.node_count());
    for (std::size_t i = 0; i < nodal.size(); ++i) nodal[i] = 1.0 + 2.0 * tri.nodes[i].x - 0.5 * tri.nodes[i].y;
    vec("fem_to_cells", multiphysics::fem_nodes_to_fvm_cells(tri, nodal, cells));
    std::vector<double> cellv(cells.cell_count());
    for (std::size_t i = 0; i < cellv.size(); ++i) cellv[i] = cells.cells()[i].center.x * cells.cells()[i].center.x;
    vec("cells_to_fem", multiphysics::fvm_cells_to_fem_nodes(cells, cellv, tri, 4));
    // Joule heating.
    {
        fem::Heat2DConfig hc; hc.conductivity = 1.5; hc.volumetric_heat_capacity = 2.0; hc.dt = 0.05;
        multiphysics::JouleHeatingCoupler2D joule(fem::make_rectangle_tri_mesh(4, 3, 0.4, 0.3), hc);
        joule.set_electrical_boundary(0, {fem::ScalarBoundaryType::dirichlet, [](fem::Node2) { return 1.0; }, {}});
        joule.set_electrical_boundary(1, {fem::ScalarBoundaryType::dirichlet, [](fem::Node2) { return 0.0; }, {}});
        joule.initialize_temperature([](fem::Node2) { return 300.0; });
        joule.set_thermal_dirichlet([](fem::Node2, double) { return 300.0; });
        auto sigma = [](fem::Node2 p) { return 2.0 + p.x; };
        joule.solve_electrical(sigma);
        joule.thermal_run(4);
        vec("joule_density", joule.joule_heating_density());
        vec("joule_temperature", joule.thermal().temperature());
    }
    // Bioheat.
    num("sar", multiphysics::sar_from_rms_electric_field(0.8, 1050.0, 40.0));
    auto fb = multiphysics::temperature_feedback({}, 312.0);
    vec("feedback", {fb.conductivity_s_per_m, fb.perfusion_per_s});
    multiphysics::VoxelTissueGrid3D vg;
    vg.nx = 3; vg.ny = 2; vg.nz = 2; vg.dx_m = vg.dy_m = vg.dz_m = 0.01;
    vg.materials = {multiphysics::VoxelTissueMaterial{}, multiphysics::VoxelTissueMaterial{1100.0, 3500.0, 0.4, 0.3, 0.0, 0.0}};
    for (std::size_t i = 0; i < 12; ++i) { vg.material_index.push_back(i % 2); vg.electric_rms_v_per_m.push_back(10.0 + 3.0 * static_cast<double>(i)); }
    auto sar = multiphysics::voxel_sar_w_per_kg(vg);
    vec("voxel_sar", sar);
    num("peak_sar_2g", multiphysics::max_mass_averaged_sar_w_per_kg(vg, sar, 0.002));
    vec("sar_projection", multiphysics::project_voxel_sar_to_pennes2d(vg, sar, 2, 2));
    for (auto boundary : {multiphysics::BioheatBoundary::periodic, multiphysics::BioheatBoundary::fixed_temperature}) {
        multiphysics::PennesBioheat2DConfig pc; pc.nx = 6; pc.ny = 5; pc.blood_perfusion_per_s = 0.005; pc.metabolic_heat_w_per_m3 = 400.0; pc.dt_s = 2.0; pc.boundary = boundary;
        multiphysics::PennesBioheat2D pennes(pc);
        std::vector<double> map(30);
        for (std::size_t i = 0; i < 30; ++i) map[i] = 0.5 * static_cast<double>(i % 7);
        pennes.set_sar(map);
        pennes.set_volumetric_heating(1000.0);
        pennes.step(3);
        vec(boundary == multiphysics::BioheatBoundary::periodic ? "pennes_periodic" : "pennes_fixed", pennes.temperature_k());
    }
    std::cout << "\n}\n";
}
